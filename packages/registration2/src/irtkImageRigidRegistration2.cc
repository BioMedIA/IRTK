/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration2.h>

#include <irtkGradientImageFilter.h>

#include <irtkHomogeneousTransformationIterator.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

void irtkImageRigidRegistration2::GuessParameter()
{
  int i;
  double xsize, ysize, zsize;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkImageRigidRegistration2::GuessParameter: Target and source image not found" << endl;
    exit(1);
  }

  // Default parameters for registration
  _NumberOfLevels     = 3;
  _NumberOfBins       = 64;

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
  _Epsilon            = 0.0001;

  // Read target pixel size
  _target->GetPixelSize(&xsize, &ysize, &zsize);

  // Default target parameters
  _TargetBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _TargetResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _TargetResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  if (_target->GetZ() > 1) {
    _TargetResolution[0][2] = GuessResolution(xsize, ysize, zsize);
  } else {
    _TargetResolution[0][2] = xsize;
  }

  for (i = 1; i < _NumberOfLevels; i++) {
    _TargetBlurring[i]      = _TargetBlurring[i-1] * 2;
    _TargetResolution[i][0] = _TargetResolution[i-1][0] * 2;
    _TargetResolution[i][1] = _TargetResolution[i-1][1] * 2;
    if (_target->GetZ() > 1) {
      _TargetResolution[i][2] = _TargetResolution[i-1][2] * 2;
    } else {
      _TargetResolution[i][2] = xsize;
    }
  }

  // Read source pixel size
  _source->GetPixelSize(&xsize, &ysize, &zsize);

  // Default source parameters
  _SourceBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _SourceResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _SourceResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  if (_source->GetZ() > 1) {
    _SourceResolution[0][2] = GuessResolution(xsize, ysize, zsize);
  } else {
    _SourceResolution[0][2] = xsize;
  }

  for (i = 1; i < _NumberOfLevels; i++) {
    _SourceBlurring[i]      = _SourceBlurring[i-1] * 2;
    _SourceResolution[i][0] = _SourceResolution[i-1][0] * 2;
    _SourceResolution[i][1] = _SourceResolution[i-1][1] * 2;
    if (_source->GetZ() > 1) {
      _SourceResolution[i][2] = _SourceResolution[i-1][2] * 2;
    } else {
      _SourceResolution[i][2] = xsize;
    }
  }

  // Remaining parameters
  for (i = 0; i < _NumberOfLevels; i++) {
    _NumberOfIterations[i] = 40;
    _MinStep[i]            = 0.01;
    _MaxStep[i]            = 1.0;
  }

  // Try to guess padding by looking at voxel values in all eight corners of the volume:
  // If all values are the same we assume that they correspond to the padding value
  _TargetPadding = MIN_GREY;
  if ((_target->Get(_target->GetX()-1, 0, 0)                                 == _target->Get(0, 0, 0)) &&
      (_target->Get(0, _target->GetY()-1, 0)                                 == _target->Get(0, 0, 0)) &&
      (_target->Get(0, 0, _target->GetZ()-1)                                 == _target->Get(0, 0, 0)) &&
      (_target->Get(_target->GetX()-1, _target->GetY()-1, 0)                 == _target->Get(0, 0, 0)) &&
      (_target->Get(0, _target->GetY()-1, _target->GetZ()-1)                 == _target->Get(0, 0, 0)) &&
      (_target->Get(_target->GetX()-1, 0, _target->GetZ()-1)                 == _target->Get(0, 0, 0)) &&
      (_target->Get(_target->GetX()-1, _target->GetY()-1, _target->GetZ()-1) == _target->Get(0, 0, 0))) {
    _TargetPadding = _target->Get(0, 0, 0);
  }
}

void irtkImageRigidRegistration2::UpdateSource()
{
  double t1, t2, u1, u2, v1, v2;
  int a, b, c, i, j, k, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

  // Pointer to reference data
  irtkRealPixel *ptr2source;
  irtkRealPixel *ptr2result;
  irtkGreyPixel  *ptr2mask;

  IRTK_START_TIMING();

  // Create iterator
  irtkHomogeneousTransformationIterator
  iterator((irtkHomogeneousTransformation *)_transformation);

  // Generate transformed tmp image
  _transformedSource = *_target;

  // Pointer to voxels in images
  ptr2result = _transformedSource.GetPointerToVoxels();
  ptr2mask   = _distanceMask.GetPointerToVoxels();

  // Calculate offsets for fast pixel access
  offset1 = 0;
  offset2 = 1;
  offset3 = this->_source->GetX();
  offset4 = this->_source->GetX()+1;
  offset5 = this->_source->GetX()*this->_source->GetY();
  offset6 = this->_source->GetX()*this->_source->GetY()+1;
  offset7 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX();
  offset8 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX()+1;

  // Initialize iterator
  iterator.Initialize(_target, _source);

  if ((_target->GetZ() == 1) && (_source->GetZ() == 1)) {

    // Loop over all voxels in the target (reference) volume
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        // Check whether reference point is valid
        if (*ptr2mask == 0) {
          // Check whether transformed point is inside source volume
          if ((iterator._x > _source_x1) && (iterator._x < _source_x2) &&
              (iterator._y > _source_y1) && (iterator._y < _source_y2)) {
            // Calculated integer coordinates
            a  = int(iterator._x);
            b  = int(iterator._y);

            // Calculated fractional coordinates
            t1 = iterator._x - a;
            u1 = iterator._y - b;
            t2 = 1 - t1;
            u2 = 1 - u1;

            // Linear interpolation in source image
            ptr2source  = _source->GetPointerToVoxels(a, b, 0);
            *ptr2result = t1 * (u2 * ptr2source[offset2] + u1 * ptr2source[offset4]) + t2 * (u2 * ptr2source[offset1] + u1 * ptr2source[offset3]);
          } else {
            *ptr2result = _SourcePadding;
          }
          iterator.NextX();
        } else {
          // Advance iterator by offset
          iterator.NextX(*ptr2mask);
          i          += (*ptr2mask) - 1;
          ptr2result += (*ptr2mask) - 1;
          ptr2mask   += (*ptr2mask) - 1;
        }
        ptr2result++;
        ptr2mask++;
      }
      iterator.NextY();
    }
  } else {
    // Loop over all voxels in the target (reference) volume
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          // Check whether reference point is valid
          if (*ptr2mask == 0) {
            // Check whether transformed point is inside source volume
            if ((iterator._x > _source_x1) && (iterator._x < _source_x2) &&
                (iterator._y > _source_y1) && (iterator._y < _source_y2) &&
                (iterator._z > _source_z1) && (iterator._z < _source_z2)) {
              // Calculated integer coordinates
              a  = int(iterator._x);
              b  = int(iterator._y);
              c  = int(iterator._z);

              // Calculated fractional coordinates
              t1 = iterator._x - a;
              u1 = iterator._y - b;
              v1 = iterator._z - c;
              t2 = 1 - t1;
              u2 = 1 - u1;
              v2 = 1 - v1;

              // Linear interpolation in source image
              ptr2source  = _source->GetPointerToVoxels(a, b, c);
              *ptr2result = (t1 * (u2 * (v2 * ptr2source[offset2] + v1 * ptr2source[offset6]) +
                                   u1 * (v2 * ptr2source[offset4] + v1 * ptr2source[offset8])) +
                             t2 * (u2 * (v2 * ptr2source[offset1] + v1 * ptr2source[offset5]) +
                                   u1 * (v2 * ptr2source[offset3] + v1 * ptr2source[offset7])));
            } else {
              *ptr2result = _SourcePadding;
            }
            iterator.NextX();
          } else {
            // Advance iterator by offset
            iterator.NextX(*ptr2mask);
            i          += (*ptr2mask) - 1;
            ptr2result += (*ptr2mask) - 1;
            ptr2mask   += (*ptr2mask) - 1;
          }
          ptr2result++;
          ptr2mask++;
        }
        iterator.NextY();
      }
      iterator.NextZ();
    }
  }

  IRTK_END_TIMING("irtkImageRigidRegistration2::UpdateSource()");
}

void irtkImageRigidRegistration2::UpdateSourceAndGradient()
{
  double t1, t2, u1, u2, v1, v2;
  int a, b, c, i, j, k, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

  // Pointer to reference data
  irtkRealPixel *ptr2source;
  irtkRealPixel *ptr2result;
  irtkRealPixel *ptr2gradX;
  irtkRealPixel *ptr2gradY;
  irtkRealPixel *ptr2gradZ;
  irtkRealPixel *ptr;
  irtkGreyPixel  *ptr2mask;

  IRTK_START_TIMING();

  // Create iterator
  irtkHomogeneousTransformationIterator
  iterator((irtkHomogeneousTransformation *)_transformation);

  // Generate transformed tmp image
  _transformedSource = *_target;

  // Pointer to voxels in images
  ptr2result = _transformedSource.GetPointerToVoxels();
  ptr2mask   = _distanceMask.GetPointerToVoxels();
  ptr2gradX  = _transformedSourceGradient.GetPointerToVoxels(0, 0, 0, 0);
  ptr2gradY  = _transformedSourceGradient.GetPointerToVoxels(0, 0, 0, 1);
  ptr2gradZ  = _transformedSourceGradient.GetPointerToVoxels(0, 0, 0, 2);

  // Calculate offsets for fast pixel access
  offset1 = 0;
  offset2 = 1;
  offset3 = this->_source->GetX();
  offset4 = this->_source->GetX()+1;
  offset5 = this->_source->GetX()*this->_source->GetY();
  offset6 = this->_source->GetX()*this->_source->GetY()+1;
  offset7 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX();
  offset8 = this->_source->GetX()*this->_source->GetY()+this->_source->GetX()+1;

  // Initialize iterator
  iterator.Initialize(_target, _source);

  if ((_target->GetZ() == 1) && (_source->GetZ() == 1)) {

    // Loop over all voxels in the target (reference) volume
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        // Check whether reference point is valid
        if (*ptr2mask == 0) {
          // Check whether transformed point is inside source volume
          if ((iterator._x > _source_x1) && (iterator._x < _source_x2) &&
              (iterator._y > _source_y1) && (iterator._y < _source_y2)) {
            // Calculated integer coordinates
            a  = int(iterator._x);
            b  = int(iterator._y);

            // Calculated fractional coordinates
            t1 = iterator._x - a;
            u1 = iterator._y - b;
            t2 = 1 - t1;
            u2 = 1 - u1;

            // Linear interpolation in source image
            ptr2source  = _source->GetPointerToVoxels(a, b, 0);
            *ptr2result = t1 * (u2 * ptr2source[offset2] + u1 * ptr2source[offset4]) + t2 * (u2 * ptr2source[offset1] + u1 * ptr2source[offset3]);

            // Linear interpolation in gradient image
            ptr = _sourceGradient.GetPointerToVoxels(a, b, 0, 0);
            *ptr2gradX = t1 * (u2 * ptr[offset2] + u1 * ptr[offset4]) + t2 * (u2 * ptr[offset1] + u1 * ptr[offset3]);
            ptr = _sourceGradient.GetPointerToVoxels(a, b, 0, 1);
            *ptr2gradY = t1 * (u2 * ptr[offset2] + u1 * ptr[offset4]) + t2 * (u2 * ptr[offset1] + u1 * ptr[offset3]);

          } else {
            *ptr2result = _SourcePadding;
            *ptr2gradX  = 0;
            *ptr2gradY  = 0;
          }
          iterator.NextX();
        } else {
          // Advance iterator by offset
          iterator.NextX(*ptr2mask);
          i          += (*ptr2mask) - 1;
          ptr2result += (*ptr2mask) - 1;
          ptr2gradX  += (*ptr2mask) - 1;
          ptr2gradY  += (*ptr2mask) - 1;
          ptr2mask   += (*ptr2mask) - 1;
        }
        ptr2result++;
        ptr2gradX++;
        ptr2gradY++;
        ptr2mask++;
      }
      iterator.NextY();
    }
  } else {
    // Loop over all voxels in the target (reference) volume
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          // Check whether reference point is valid
          if (*ptr2mask == 0) {
            // Check whether transformed point is inside source volume
            if ((iterator._x > _source_x1) && (iterator._x < _source_x2) &&
                (iterator._y > _source_y1) && (iterator._y < _source_y2) &&
                (iterator._z > _source_z1) && (iterator._z < _source_z2)) {
              // Calculated integer coordinates
              a  = int(iterator._x);
              b  = int(iterator._y);
              c  = int(iterator._z);

              // Calculated fractional coordinates
              t1 = iterator._x - a;
              u1 = iterator._y - b;
              v1 = iterator._z - c;
              t2 = 1 - t1;
              u2 = 1 - u1;
              v2 = 1 - v1;

              // Linear interpolation in source image
              ptr2source  = _source->GetPointerToVoxels(a, b, c);
              *ptr2result = (t1 * (u2 * (v2 * ptr2source[offset2] + v1 * ptr2source[offset6]) +
                                   u1 * (v2 * ptr2source[offset4] + v1 * ptr2source[offset8])) +
                             t2 * (u2 * (v2 * ptr2source[offset1] + v1 * ptr2source[offset5]) +
                                   u1 * (v2 * ptr2source[offset3] + v1 * ptr2source[offset7])));


              // Linear interpolation in gradient image
              ptr = _sourceGradient.GetPointerToVoxels(a, b, c, 0);
              *ptr2gradX = (t1 * (u2 * (v2 * ptr[offset2] + v1 * ptr[offset6]) +
                                  u1 * (v2 * ptr[offset4] + v1 * ptr[offset8])) +
                            t2 * (u2 * (v2 * ptr[offset1] + v1 * ptr[offset5]) +
                                  u1 * (v2 * ptr[offset3] + v1 * ptr[offset7])));
              ptr = _sourceGradient.GetPointerToVoxels(a, b, c, 1);
              *ptr2gradY = (t1 * (u2 * (v2 * ptr[offset2] + v1 * ptr[offset6]) +
                                  u1 * (v2 * ptr[offset4] + v1 * ptr[offset8])) +
                            t2 * (u2 * (v2 * ptr[offset1] + v1 * ptr[offset5]) +
                                  u1 * (v2 * ptr[offset3] + v1 * ptr[offset7])));
              ptr = _sourceGradient.GetPointerToVoxels(a, b, c, 2);
              *ptr2gradZ = (t1 * (u2 * (v2 * ptr[offset2] + v1 * ptr[offset6]) +
                                  u1 * (v2 * ptr[offset4] + v1 * ptr[offset8])) +
                            t2 * (u2 * (v2 * ptr[offset1] + v1 * ptr[offset5]) +
                                  u1 * (v2 * ptr[offset3] + v1 * ptr[offset7])));

            } else {
              *ptr2result = _SourcePadding;
              *ptr2gradX  = 0;
              *ptr2gradY  = 0;
              *ptr2gradZ  = 0;
            }
            iterator.NextX();
          } else {
            // Advance iterator by offset
            iterator.NextX(*ptr2mask);
            i          += (*ptr2mask) - 1;
            ptr2result += (*ptr2mask) - 1;
            ptr2gradX  += (*ptr2mask) - 1;
            ptr2gradY  += (*ptr2mask) - 1;
            ptr2gradZ  += (*ptr2mask) - 1;
            ptr2mask   += (*ptr2mask) - 1;
          }
          ptr2result++;
          ptr2gradX++;
          ptr2gradY++;
          ptr2gradZ++;
          ptr2mask++;
        }
        iterator.NextY();
      }
      iterator.NextZ();
    }
  }

  IRTK_END_TIMING("irtkImageRigidRegistration2::UpdateSourceAndGradient()");
}

double irtkImageRigidRegistration2::EvaluateGradient(double *gradient)
{
  int i, j, k, l;
  double x, y, z, jac[3], norm, max_length;
  static double *g = NULL, *h = NULL, gg, dgg, gamma;

  // Pointer to reference data
  irtkGreyPixel  *ptr2mask;
  irtkRealPixel *ptr2gradX, *ptr2gradY, *ptr2gradZ;

  // Compute gradient with respect to displacements
  this->irtkImageRegistration2::EvaluateGradient(gradient);

  IRTK_START_TIMING();

  // Convert this gradient into gradient with respect to parameters
  for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
    gradient[i] = 0;
  }

  // Pointer to voxels in images
  ptr2mask   = _distanceMask.GetPointerToVoxels();
  ptr2gradX  = _similarityGradient.GetPointerToVoxels(0, 0, 0, 0);
  ptr2gradY  = _similarityGradient.GetPointerToVoxels(0, 0, 0, 1);
  ptr2gradZ  = _similarityGradient.GetPointerToVoxels(0, 0, 0, 2);

  // Loop over images
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        // Check whether reference point is valid
        if (*ptr2mask == 0) {
          x = i;
          y = j;
          z = k;

          // Convert voxel coordinate to world coordinate
          _target->ImageToWorld(x, y, z);

          // Convert this gradient into gradient with respect to parameters
          for (l = 0; l < _transformation->NumberOfDOFs(); l++) {

            // Compute derivatives with respect to DOF
            _transformation->JacobianDOFs(jac, l, x, y, z);
            gradient[l] += jac[0] * *ptr2gradX + jac[1] * *ptr2gradY + jac[2] * *ptr2gradZ;
          }
        }
        ptr2gradX++;
        ptr2gradY++;
        ptr2gradZ++;
        ptr2mask++;
      }
    }
  }

  // Update gradient
  if (_CurrentIteration == 0) {
    // First iteration, so let's initialize
    if (g != NULL) delete []g;
    g = new double [_transformation->NumberOfDOFs()];
    if (h != NULL) delete []h;
    h = new double [_transformation->NumberOfDOFs()];
    for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
      g[i] = -gradient[i];
      h[i] = g[i];
    }
  } else {
    // Update gradient direction to be conjugate
    gg = 0;
    dgg = 0;
    for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
      gg  += g[i]*h[i];
      dgg += (gradient[i]+g[i])*gradient[i];
    }
    gamma = dgg/gg;
    for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
      g[i] = -gradient[i];
      h[i] = g[i] + gamma*h[i];
      gradient[i] = -h[i];
    }
  }

  // Calculate maximum of gradient vector
  max_length = 0;
  for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
    norm = sqrt(gradient[i] * gradient[i]);
    if (norm > max_length) max_length = norm;
  }

  // Deal with active and passive DOFs
  for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
    if (_transformation->irtkTransformation::GetStatus(i) == _Passive) {
      gradient[i] = 0;
    }
  }

  IRTK_END_TIMING("irtkImageRigidRegistration2::EvaluateGradient()");

  return max_length;
}

