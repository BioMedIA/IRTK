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

void irtkMultipleImageRigidRegistration2::GuessParameter()
{
  int i;

  for(i = 0; i < _numberOfImages; i++){
      if ((_target[i] == NULL) || (_source[i] == NULL)) {
          cerr << "irtkMultipleImageRigidRegistration2::GuessParameter: Target and source image "<<i<<" not found" << endl;
          exit(1);
      }
  }

  // Default parameters for registration
  _NumberOfLevels     = 3;
  _NumberOfBins       = 64;

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
  _Epsilon            = 0.0001;
  _Lregu              = 0.0001;

  // Default target parameters
  _TargetBlurring[0]      = 0.5;

  for (i = 1; i < _NumberOfLevels; i++) {
    _TargetBlurring[i]      = _TargetBlurring[i-1] * 2;
  }

  // Default source parameters
  _SourceBlurring[0]      = 0.5;

  for (i = 1; i < _NumberOfLevels; i++) {
    _SourceBlurring[i]      = _SourceBlurring[i-1] * 2;
  }

  // Remaining parameters
  for (i = 0; i < _NumberOfLevels; i++) {
    _NumberOfIterations[i] = 40;
    _MinStep[i]            = 0.01;
    _MaxStep[i]            = 2.0;
  }

  // Try to guess padding by looking at voxel values in all eight corners of the volume:
  // If all values are the same we assume that they correspond to the padding value
  _TargetPadding = MIN_GREY;
}

void irtkMultipleImageRigidRegistration2::UpdateSource()
{
  double t1, t2, u1, u2, v1, v2;
  int a, b, c, i, j, k, n, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

  // Pointer to reference data
  short *ptr2target;
  short *ptr2source;
  double *ptr2result;

  // Start timing
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  for(n = 0; n < _numberOfImages; n++){

      // Create iterator
      irtkHomogeneousTransformationIterator
          iterator((irtkHomogeneousTransformation *)_transformation);

      // Generate transformed tmp image
      _transformedSource[n] = *_target[n];

      // Pointer to voxels in images
      ptr2target = _target[n]->GetPointerToVoxels();
      ptr2result = _transformedSource[n].GetPointerToVoxels();

      // Calculate offsets for fast pixel access
      offset1 = 0;
      offset2 = 1;
      offset3 = this->_source[n]->GetX();
      offset4 = this->_source[n]->GetX()+1;
      offset5 = this->_source[n]->GetX()*this->_source[n]->GetY();
      offset6 = this->_source[n]->GetX()*this->_source[n]->GetY()+1;
      offset7 = this->_source[n]->GetX()*this->_source[n]->GetY()+this->_source[n]->GetX();
      offset8 = this->_source[n]->GetX()*this->_source[n]->GetY()+this->_source[n]->GetX()+1;

      // Initialize iterator
      iterator.Initialize(_target[n], _source[n]);

      if ((_target[n]->GetZ() == 1) && (_source[n]->GetZ() == 1)) {

          // Loop over all voxels in the target (reference) volume
          for (j = 0; j < _target[n]->GetY(); j++) {
              for (i = 0; i < _target[n]->GetX(); i++) {
                  // Check whether reference point is valid
                  if (*ptr2target >= 0) {
                      // Check whether transformed point is inside source volume
                      if ((iterator._x > _source_x1[n]) && (iterator._x < _source_x2[n]) &&
                          (iterator._y > _source_y1[n]) && (iterator._y < _source_y2[n])) {
                              // Calculated integer coordinates
                              a  = int(iterator._x);
                              b  = int(iterator._y);

                              // Calculated fractional coordinates
                              t1 = iterator._x - a;
                              u1 = iterator._y - b;
                              t2 = 1 - t1;
                              u2 = 1 - u1;

                              // Linear interpolation in source image
                              ptr2source  = _source[n]->GetPointerToVoxels(a, b, 0);
                              *ptr2result = t1 * (u2 * ptr2source[offset2] + u1 * ptr2source[offset4]) + t2 * (u2 * ptr2source[offset1] + u1 * ptr2source[offset3]);
                      } else {
                          *ptr2result = -1;
                      }
                      iterator.NextX();
                  } else {
                      // Advance iterator by offset
                      iterator.NextX(*ptr2target * -1);
                      i          -= (*ptr2target) + 1;
                      ptr2result -= (*ptr2target) + 1;
                      ptr2target -= (*ptr2target) + 1;
                  }
                  ptr2target++;
                  ptr2result++;
              }
              iterator.NextY();
          }
      } else {
          // Loop over all voxels in the target (reference) volume
          for (k = 0; k < _target[n]->GetZ(); k++) {
              for (j = 0; j < _target[n]->GetY(); j++) {
                  for (i = 0; i < _target[n]->GetX(); i++) {
                      // Check whether reference point is valid
                      if (*ptr2target >= 0) {
                          // Check whether transformed point is inside source volume
                          if ((iterator._x > _source_x1[n]) && (iterator._x < _source_x2[n]) &&
                              (iterator._y > _source_y1[n]) && (iterator._y < _source_y2[n]) &&
                              (iterator._z > _source_z1[n]) && (iterator._z < _source_z2[n])) {
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
                                  ptr2source  = _source[n]->GetPointerToVoxels(a, b, c);
                                  *ptr2result = (t1 * (u2 * (v2 * ptr2source[offset2] + v1 * ptr2source[offset6]) +
                                      u1 * (v2 * ptr2source[offset4] + v1 * ptr2source[offset8])) +
                                      t2 * (u2 * (v2 * ptr2source[offset1] + v1 * ptr2source[offset5]) +
                                      u1 * (v2 * ptr2source[offset3] + v1 * ptr2source[offset7])));
                          } else {
                              *ptr2result = -1;
                          }
                          iterator.NextX();
                      } else {
                          // Advance iterator by offset
                          iterator.NextX(*ptr2target * -1);
                          i          -= (*ptr2target) + 1;
                          ptr2result -= (*ptr2target) + 1;
                          ptr2target -= (*ptr2target) + 1;
                      }
                      ptr2target++;
                      ptr2result++;
                  }
                  iterator.NextY();
              }
              iterator.NextZ();
          }
      }
  }

  // Stop timing
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //cout << "CPU time for irtkMultipleImageRigidRegistration2::Update() = " << cpu_time_used << endl;
}

void irtkMultipleImageRigidRegistration2::UpdateSourceAndGradient()
{
  double t1, t2, u1, u2, v1, v2;
  int a, b, c, i, j, k, n, offset1, offset2, offset3, offset4, offset5, offset6, offset7, offset8;

  // Pointer to reference data
  short *ptr2target;
  short *ptr2source;
  double *ptr2result;
  double *ptr2gradX;
  double *ptr2gradY;
  double *ptr2gradZ;
  double *ptr;

  // Start timing
  clock_t start, end;
  double cpu_time_used;
  start = clock();
  for(n = 0; n < _numberOfImages; n++){
      // Create iterator
      irtkHomogeneousTransformationIterator
          iterator((irtkHomogeneousTransformation *)_transformation);

      // Generate transformed tmp image
      _transformedSource[n] = *_target[n];

      // Pointer to voxels in images
      ptr2target = _target[n]->GetPointerToVoxels();
      ptr2result = _transformedSource[n].GetPointerToVoxels();
      ptr2gradX  = _transformedSourceGradient[n].GetPointerToVoxels(0, 0, 0, 0);
      ptr2gradY  = _transformedSourceGradient[n].GetPointerToVoxels(0, 0, 0, 1);
      ptr2gradZ  = _transformedSourceGradient[n].GetPointerToVoxels(0, 0, 0, 2);

      // Calculate offsets for fast pixel access
      offset1 = 0;
      offset2 = 1;
      offset3 = this->_source[n]->GetX();
      offset4 = this->_source[n]->GetX()+1;
      offset5 = this->_source[n]->GetX()*this->_source[n]->GetY();
      offset6 = this->_source[n]->GetX()*this->_source[n]->GetY()+1;
      offset7 = this->_source[n]->GetX()*this->_source[n]->GetY()+this->_source[n]->GetX();
      offset8 = this->_source[n]->GetX()*this->_source[n]->GetY()+this->_source[n]->GetX()+1;

      // Initialize iterator
      iterator.Initialize(_target[n], _source[n]);

      if ((_target[n]->GetZ() == 1) && (_source[n]->GetZ() == 1)) {

          // Loop over all voxels in the target (reference) volume
          for (j = 0; j < _target[n]->GetY(); j++) {
              for (i = 0; i < _target[n]->GetX(); i++) {
                  // Check whether reference point is valid
                  if (*ptr2target >= 0) {
                      // Check whether transformed point is inside source volume
                      if ((iterator._x > _source_x1[n]) && (iterator._x < _source_x2[n]) &&
                          (iterator._y > _source_y1[n]) && (iterator._y < _source_y2[n])) {
                              // Calculated integer coordinates
                              a  = int(iterator._x);
                              b  = int(iterator._y);

                              // Calculated fractional coordinates
                              t1 = iterator._x - a;
                              u1 = iterator._y - b;
                              t2 = 1 - t1;
                              u2 = 1 - u1;

                              // Linear interpolation in source image
                              ptr2source  = _source[n]->GetPointerToVoxels(a, b, 0);
                              *ptr2result = t1 * (u2 * ptr2source[offset2] + u1 * ptr2source[offset4]) + t2 * (u2 * ptr2source[offset1] + u1 * ptr2source[offset3]);

                              // Linear interpolation in gradient image
                              ptr = _sourceGradient[n].GetPointerToVoxels(a, b, 0, 0);
                              *ptr2gradX = t1 * (u2 * ptr[offset2] + u1 * ptr[offset4]) + t2 * (u2 * ptr[offset1] + u1 * ptr[offset3]);
                              ptr = _sourceGradient[n].GetPointerToVoxels(a, b, 0, 1);
                              *ptr2gradY = t1 * (u2 * ptr[offset2] + u1 * ptr[offset4]) + t2 * (u2 * ptr[offset1] + u1 * ptr[offset3]);

                      } else {
                          *ptr2result = -1;
                          *ptr2gradX  = 0;
                          *ptr2gradY  = 0;
                      }
                      iterator.NextX();
                  } else {
                      // Advance iterator by offset
                      iterator.NextX(*ptr2target * -1);
                      i          -= (*ptr2target) + 1;
                      ptr2result -= (*ptr2target) + 1;
                      ptr2gradX  -= (*ptr2target) + 1;
                      ptr2gradY  -= (*ptr2target) + 1;
                      ptr2target -= (*ptr2target) + 1;
                  }
                  ptr2target++;
                  ptr2result++;
                  ptr2gradX++;
                  ptr2gradY++;
              }
              iterator.NextY();
          }
      } else {
          // Loop over all voxels in the target (reference) volume
          for (k = 0; k < _target[n]->GetZ(); k++) {
              for (j = 0; j < _target[n]->GetY(); j++) {
                  for (i = 0; i < _target[n]->GetX(); i++) {
                      // Check whether reference point is valid
                      if (*ptr2target >= 0) {
                          // Check whether transformed point is inside source volume
                          if ((iterator._x > _source_x1[n]) && (iterator._x < _source_x2[n]) &&
                              (iterator._y > _source_y1[n]) && (iterator._y < _source_y2[n]) &&
                              (iterator._z > _source_z1[n]) && (iterator._z < _source_z2[n])) {
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
                                  ptr2source  = _source[n]->GetPointerToVoxels(a, b, c);
                                  *ptr2result = (t1 * (u2 * (v2 * ptr2source[offset2] + v1 * ptr2source[offset6]) +
                                      u1 * (v2 * ptr2source[offset4] + v1 * ptr2source[offset8])) +
                                      t2 * (u2 * (v2 * ptr2source[offset1] + v1 * ptr2source[offset5]) +
                                      u1 * (v2 * ptr2source[offset3] + v1 * ptr2source[offset7])));


                                  // Linear interpolation in gradient image
                                  ptr = _sourceGradient[n].GetPointerToVoxels(a, b, c, 0);
                                  *ptr2gradX = (t1 * (u2 * (v2 * ptr[offset2] + v1 * ptr[offset6]) +
                                      u1 * (v2 * ptr[offset4] + v1 * ptr[offset8])) +
                                      t2 * (u2 * (v2 * ptr[offset1] + v1 * ptr[offset5]) +
                                      u1 * (v2 * ptr[offset3] + v1 * ptr[offset7])));
                                  ptr = _sourceGradient[n].GetPointerToVoxels(a, b, c, 1);
                                  *ptr2gradY = (t1 * (u2 * (v2 * ptr[offset2] + v1 * ptr[offset6]) +
                                      u1 * (v2 * ptr[offset4] + v1 * ptr[offset8])) +
                                      t2 * (u2 * (v2 * ptr[offset1] + v1 * ptr[offset5]) +
                                      u1 * (v2 * ptr[offset3] + v1 * ptr[offset7])));
                                  ptr = _sourceGradient[n].GetPointerToVoxels(a, b, c, 2);
                                  *ptr2gradZ = (t1 * (u2 * (v2 * ptr[offset2] + v1 * ptr[offset6]) +
                                      u1 * (v2 * ptr[offset4] + v1 * ptr[offset8])) +
                                      t2 * (u2 * (v2 * ptr[offset1] + v1 * ptr[offset5]) +
                                      u1 * (v2 * ptr[offset3] + v1 * ptr[offset7])));

                          } else {
                              *ptr2result = -1;
                              *ptr2gradX  = 0;
                              *ptr2gradY  = 0;
                              *ptr2gradZ  = 0;
                          }
                          iterator.NextX();
                      } else {
                          // Advance iterator by offset
                          iterator.NextX(*ptr2target * -1);
                          i          -= (*ptr2target) + 1;
                          ptr2result -= (*ptr2target) + 1;
                          ptr2gradX  -= (*ptr2target) + 1;
                          ptr2gradY  -= (*ptr2target) + 1;
                          ptr2gradZ  -= (*ptr2target) + 1;
                          ptr2target -= (*ptr2target) + 1;
                      }
                      ptr2target++;
                      ptr2result++;
                      ptr2gradX++;
                      ptr2gradY++;
                      ptr2gradZ++;
                  }
                  iterator.NextY();
              }
              iterator.NextZ();
          }
      }
  }

  // Stop timing
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //cout << "CPU time for irtkMultipleImageRigidRegistration2::Update() = " << cpu_time_used << endl;
}

double irtkMultipleImageRigidRegistration2::EvaluateGradient(double *gradient)
{
  int i, j, k, l, n;
  double x, y, z, jac[3], norm, max_length;
  static double *g = NULL, *h = NULL, gg, dgg, gamma;

  // Pointer to reference data
  short *ptr2target;
  double *ptr2result, *ptr2gradX, *ptr2gradY, *ptr2gradZ;

  // Compute gradient with respect to displacements
  this->irtkMultipleImageRegistration2::EvaluateGradient(gradient);

  // Start timing
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  // Convert this gradient into gradient with respect to parameters
  for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
    gradient[i] = 0;
  }

  for(n = 0; n < _numberOfImages; n++){

      // Pointer to voxels in images
      ptr2target = _target[n]->GetPointerToVoxels();
      ptr2result = _transformedSource[n].GetPointerToVoxels();
      ptr2gradX  = _similarityGradient[n].GetPointerToVoxels(0, 0, 0, 0);
      ptr2gradY  = _similarityGradient[n].GetPointerToVoxels(0, 0, 0, 1);
      ptr2gradZ  = _similarityGradient[n].GetPointerToVoxels(0, 0, 0, 2);

      // Loop over images
      for (k = 0; k < _target[n]->GetZ(); k++) {
          for (j = 0; j < _target[n]->GetY(); j++) {
              for (i = 0; i < _target[n]->GetX(); i++) {
                  // Check whether reference point is valid
                  if ((*ptr2target >= 0) && (*ptr2result >= 0)) {
                      x = i;
                      y = j;
                      z = k;

                      // Convert voxel coordinate to world coordinate
                      _target[n]->ImageToWorld(x, y, z);

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
                  ptr2target++;
                  ptr2result++;
              }
          }
      }
  }

#ifdef HAS_VTK
  if (_ptarget != NULL && _psource != NULL && _Lregu > 0){
      double p[3],g[3];
      double norm;
      norm = _totalVoxels/_transformation->NumberOfDOFs()*_Lregu;
      for(n = 0; n < _pgradient->GetNumberOfPoints(); n++){

          _ptarget->GetPoints()->GetPoint(n,p);

          _pgradient->GetPoints()->GetPoint(n,g);

          // Convert this gradient into gradient with respect to parameters
          for (l = 0; l < _transformation->NumberOfDOFs(); l++) {

              // Compute derivatives with respect to DOF
              _transformation->JacobianDOFs(jac, l, p[0], p[1], p[2]);
              gradient[l] += (jac[0] * g[0] + jac[1] * g[1] + jac[2] * g[2]) * norm;
          }
      }
  }
#endif

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

  // Stop timing
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //cout << "CPU time for irtkMultipleImageRigidRegistration2::EvaluateGradient() = " << cpu_time_used << endl;

  return max_length;
}

