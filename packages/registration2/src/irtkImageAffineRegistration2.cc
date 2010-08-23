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

double irtkImageAffineRegistration2::EvaluateGradient(double *gradient)
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
  cout << "CPU time for irtkImageAffineRegistration2::EvaluateGradient() = " << cpu_time_used << endl;

  return norm;
}

