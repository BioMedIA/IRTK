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
  _Epsilon            = 0.0001;

  // Read target pixel size
  _target->GetPixelSize(&xsize, &ysize, &zsize);

  // Default target parameters
  _TargetBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _TargetResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _TargetResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  _TargetResolution[0][2] = GuessResolution(xsize, ysize, zsize);

  for (i = 1; i < _NumberOfLevels; i++) {
    _TargetBlurring[i]      = _TargetBlurring[i-1] * 2;
    _TargetResolution[i][0] = _TargetResolution[i-1][0] * 2;
    _TargetResolution[i][1] = _TargetResolution[i-1][1] * 2;
    _TargetResolution[i][2] = _TargetResolution[i-1][2] * 2;
  }

  // Read source pixel size
  _source->GetPixelSize(&xsize, &ysize, &zsize);

  // Default source parameters
  _SourceBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _SourceResolution[0][0] = GuessResolution(xsize, ysize, zsize);
  _SourceResolution[0][1] = GuessResolution(xsize, ysize, zsize);
  _SourceResolution[0][2] = GuessResolution(xsize, ysize, zsize);

  for (i = 1; i < _NumberOfLevels; i++) {
    _SourceBlurring[i]      = _SourceBlurring[i-1] * 2;
    _SourceResolution[i][0] = _SourceResolution[i-1][0] * 2;
    _SourceResolution[i][1] = _SourceResolution[i-1][1] * 2;
    _SourceResolution[i][2] = _SourceResolution[i-1][2] * 2;
  }

  // Remaining parameters
  for (i = 0; i < _NumberOfLevels; i++) {
    _NumberOfIterations[i] = 20;
    _NumberOfSteps[i]      = 4;
    _LengthOfSteps[i]      = 2 * pow(2.0, i);
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

void irtkImageRigidRegistration2::Initialize()
{
  // Call base class
  this->irtkImageRegistration2::Initialize();

  // Invert rigid transformation (to be backwards compatible)
  ((irtkRigidTransformation *)_transformation)->Invert();
  ((irtkRigidTransformation *)_transformation)->UpdateParameter();
}

void irtkImageRigidRegistration2::Finalize()
{
  // Call base class
  this->irtkImageRegistration2::Finalize();

  // Invert rigid transformation (to be backwards compatible)
  ((irtkRigidTransformation *)_transformation)->Invert();
  ((irtkRigidTransformation *)_transformation)->UpdateParameter();
}

void irtkImageRigidRegistration2::Update()
{
  int i, j, k, t;

  // Pointer to reference data
  short *ptr2target;
  double *ptr2result;

  // Start timing
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  // Invert transformation
  ((irtkRigidTransformation *)_transformation)->Invert();

  // Create iterator
  irtkHomogeneousTransformationIterator
  iterator((irtkHomogeneousTransformation *)_transformation);

  // Generate transformed tmp image
  _transformedSource = *_target;

  // Pointer to voxels in images
  ptr2target = _target->GetPointerToVoxels();
  ptr2result = _transformedSource.GetPointerToVoxels();

  for (t = 0; t < _target->GetT(); t++) {

    // Initialize iterator
    iterator.Initialize(_target, _source);

    // Loop over all voxels in the target (reference) volume
    for (k = 0; k < _target->GetZ(); k++) {
      for (j = 0; j < _target->GetY(); j++) {
        for (i = 0; i < _target->GetX(); i++) {
          // Check whether reference point is valid
          if (*ptr2target >= 0) {
            // Check whether transformed point is inside source volume
            if ((iterator._x > _source_x1) && (iterator._x < _source_x2) &&
                (iterator._y > _source_y1) && (iterator._y < _source_y2) &&
                (iterator._z > _source_z1) && (iterator._z < _source_z2)) {
              // Add sample to metric
              *ptr2result = _interpolator->EvaluateInside(iterator._x, iterator._y, iterator._z, t);
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

  // Invert transformation
  ((irtkRigidTransformation *)_transformation)->Invert();

  // Compute gradient of source image
  irtkGradientImageFilter<double> gradient(irtkGradientImageFilter<double>::GRADIENT_VECTOR);
  gradient.SetInput (&_transformedSource);
  gradient.SetOutput(&_transformedSourceGradient);
  gradient.Run();

  // Stop timing
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //cout << "CPU time for irtkImageRigidRegistration2::Update() = " << cpu_time_used << endl;
}

double irtkImageRigidRegistration2::EvaluateGradient(double *gradient)
{
  int i, j, k, l;
  double x, y, z, jac [3], norm;

  // Pointer to reference data
  short *ptr2target;
  double *ptr2result, *ptr2gradX, *ptr2gradY, *ptr2gradZ;

  // Compute gradient with respect to displacements
  this->irtkImageRegistration2::EvaluateGradient(gradient);

  // Start timing
  clock_t start, end;
  double cpu_time_used;
  start = clock();

  // Convert this gradient into gradient with respect to parameters
  for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
  	gradient[i] = 0;
  }

  // Pointer to voxels in images
  ptr2target = _target->GetPointerToVoxels();
  ptr2result = _transformedSource.GetPointerToVoxels();
  ptr2gradX  = _similarityGradient.GetPointerToVoxels(0, 0, 0, 0);
  ptr2gradY  = _similarityGradient.GetPointerToVoxels(0, 0, 0, 1);
  ptr2gradZ  = _similarityGradient.GetPointerToVoxels(0, 0, 0, 2);

  // Loop over images
  for (k = 0; k < _target->GetZ(); k++) {
    for (j = 0; j < _target->GetY(); j++) {
      for (i = 0; i < _target->GetX(); i++) {
        // Check whether reference point is valid
        if ((*ptr2target >= 0) && (*ptr2result >= 0)) {
          x = i;
          y = j;
          z = k;

          // Convert voxel coordinate to world coordinate
          _target->ImageToWorld(x, y, z);

          // Convert this gradient into gradient with respect to parameters (note that the sign of the gradient is changed here)
          for (l = 0; l < _transformation->NumberOfDOFs(); l++) {

          	// Compute derivatives with respect to DOF
            _transformation->JacobianDOFs(jac, l, x, y, z);
            gradient[l] -= jac[0] * *ptr2gradX + jac[1] * *ptr2gradY + jac[2] * *ptr2gradZ;
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

  // Calculate norm of vector
  norm = 0;
  for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
    norm += gradient[i] * gradient[i];
  }

  // Normalize vector
  norm = sqrt(norm);
  if (norm > 0) {
    for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
    	gradient[i] /= norm;
    }
  } else {
    for (i = 0; i < _transformation->NumberOfDOFs(); i++) {
    	gradient[i] = 0;
    }
  }

  // Stop timing
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  //cout << "CPU time for irtkImageRigidRegistration2::EvaluateGradient() = " << cpu_time_used << endl;

  return norm;
}

