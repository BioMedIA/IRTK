/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

#include <irtkHomogeneousTransformationIterator.h>

#include <irtkMultiThreadedImageRigidRegistration.h>

void irtkImageRigidRegistration::GuessParameter()
{
  int i;
  double xsize, ysize, zsize;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkImageRigidRegistration::GuessParameter: Target and source image not found" << endl;
    exit(1);
  }

  // Default parameters for registration
  _NumberOfLevels     = 3;
  _NumberOfBins       = 64;

  // Default parameters for optimization
  _SimilarityMeasure  = NMI;
  _OptimizationMethod = GradientDescent;
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

void irtkImageRigidRegistration::Initialize()
{
  // Call base class
  this->irtkImageRegistration::Initialize();

  // Invert rigid transformation (to be backwards compatible)
  //((irtkRigidTransformation *)_transformation)->Invert();
  //((irtkRigidTransformation *)_transformation)->UpdateParameter();
}

void irtkImageRigidRegistration::Finalize()
{
  // Call base class
  this->irtkImageRegistration::Finalize();

  // Invert rigid transformation (to be backwards compatible)
  //((irtkRigidTransformation *)_transformation)->Invert();
  //((irtkRigidTransformation *)_transformation)->UpdateParameter();
}

double irtkImageRigidRegistration::Evaluate()
{

#ifndef HAS_TBB
  int i, j, k, t;
#endif

  // Pointer to reference data
  irtkGreyPixel *ptr2target;

  // Print debugging information
  this->Debug("irtkImageRigidRegistration::Evaluate");

  // Invert transformation
  //((irtkRigidTransformation *)_transformation)->Invert();

  // Create iterator
  irtkHomogeneousTransformationIterator
  iterator((irtkHomogeneousTransformation *)_transformation);

  // Initialize metric
  _metric->Reset();

  // Pointer to voxels in target image
  ptr2target = _target->GetPointerToVoxels();

#ifdef HAS_TBB
  irtkMultiThreadedImageRigidRegistrationEvaluate evaluate(this);
  parallel_reduce(blocked_range<int>(0, _target->GetZ(), 20), evaluate);
#else

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
              _metric->Add(*ptr2target, round(_interpolator->EvaluateInside(iterator._x, iterator._y, iterator._z, t)));
            }
            iterator.NextX();
          } else {
            // Advance iterator by offset
            iterator.NextX(*ptr2target * -1);
            i          -= (*ptr2target) + 1;
            ptr2target -= (*ptr2target) + 1;
          }
          ptr2target++;
        }
        iterator.NextY();
      }
      iterator.NextZ();
    }
  }

#endif


  // Invert transformation
  //((irtkRigidTransformation *)_transformation)->Invert();

  // Evaluate similarity measure
  return _metric->Evaluate();
}
