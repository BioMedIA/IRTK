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

// Used as temporary memory for transformed intensities
extern irtkGreyImage *_tmpImage;

// Used as temporary memory for metric
extern irtkSimilarityMetric *_tmpMetricA, *_tmpMetricB;

// Used as lookup table for transformed coordinates up to level n-1.  This
// lookup table needs to be calculated only once for each image resolution
// level.
extern float *_mffdLookupTable;

// Used as lookup table for transformed coordinates including level n. This
// lookup table needs to be calculated each time a control point has been
// modified.
extern float *_affdLookupTable;

// Used as lookup table for the contribution of each control point. This
// lookup table needs to be calculated only once.
extern float *_localLookupTable;

void irtkImageFreeFormRegistration2D::GuessParameter()
{
  int i;
  double xsize, ysize, zsize, spacing;

  if ((_target == NULL) || (_source == NULL)) {
    cerr << "irtkImageFreeFormRegistration2D::GuessParameter: Target and source image not found" << endl;
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

  // Use xsize as spacing
  spacing = xsize;

  // Default target parameters
  _TargetBlurring[0]      = GuessResolution(xsize, ysize) / 2.0;
  _TargetResolution[0][0] = GuessResolution(xsize, ysize);
  _TargetResolution[0][1] = GuessResolution(xsize, ysize);
  _TargetResolution[0][2] = zsize;

  for (i = 1; i < _NumberOfLevels; i++) {
    _TargetBlurring[i]      = _TargetBlurring[i-1] * 2;
    _TargetResolution[i][0] = _TargetResolution[i-1][0] * 2;
    _TargetResolution[i][1] = _TargetResolution[i-1][1] * 2;
    _TargetResolution[i][2] = zsize;
  }

  // Read source pixel size
  _source->GetPixelSize(&xsize, &ysize, &zsize);

  // Default source parameters
  _SourceBlurring[0]      = GuessResolution(xsize, ysize, zsize) / 2.0;
  _SourceResolution[0][0] = GuessResolution(xsize, ysize);
  _SourceResolution[0][1] = GuessResolution(xsize, ysize);
  _SourceResolution[0][2] = zsize;

  for (i = 1; i < _NumberOfLevels; i++) {
    _SourceBlurring[i]      = _SourceBlurring[i-1] * 2;
    _SourceResolution[i][0] = _SourceResolution[i-1][0] * 2;
    _SourceResolution[i][1] = _SourceResolution[i-1][1] * 2;
    _SourceResolution[i][2] = zsize;
  }

  // Default parameters for non-rigid registration
  _Lambda1            = 0;
  _Lambda2            = 0;
  _Lambda3            = 0;
  _DX                 =_target->GetX() * spacing / 10.0;
  _DY                 =_target->GetX() * spacing / 10.0;
  _DZ                 = 1;
  _Subdivision        = true;

  // Remaining parameters
  for (i = 0; i < _NumberOfLevels; i++) {
    _NumberOfIterations[i] = 10;
    _NumberOfSteps[i]      = 4;
    _LengthOfSteps[i]      = _DX / 8.0 * pow(2.0, i);
  }

  // Try to guess padding by looking at voxel values in all eight corners of the volume:
  // If all values are the same we assume that they correspond to the padding value
  _TargetPadding = MIN_GREY;
  if ((_target->Get(_target->GetX()-1, 0, 0)                 == _target->Get(0, 0, 0)) &&
      (_target->Get(0, _target->GetY()-1, 0)                 == _target->Get(0, 0, 0)) &&
      (_target->Get(_target->GetX()-1, _target->GetY()-1, 0) == _target->Get(0, 0, 0))) {
    _TargetPadding = _target->Get(0, 0, 0);
  }
}

double irtkImageFreeFormRegistration2D::Evaluate()
{
  // Image coordinates
  int i, j;
  // World coordinates
  double x, y, z;
  // Pointer to reference data
  irtkGreyPixel *ptr2target;
  irtkGreyPixel *ptr2tmp;
  float *ptr;

  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration2D::Evaluate");

  // Initialize histogram
  _metric->Reset();

  // Loop over all voxels in the target (reference) volume
  ptr2target = _target->GetPointerToVoxels();
  ptr2tmp    = _tmpImage->GetPointerToVoxels();
  ptr        = _mffdLookupTable;;
  for (j = 0; j < _target->GetY(); j++) {
    for (i = 0; i < _target->GetX(); i++) {
      // Check whether reference point is valid
      if (*ptr2target >= 0) {
        x = i;
        y = j;
        z = 0;
        _target->ImageToWorld(x, y, z);
        _affd->LocalDisplacement(x, y, z);
        x += ptr[0];
        y += ptr[1];
        _source->WorldToImage(x, y, z);
        // Check whether transformed point is inside volume
        if ((x > _source_x1) && (x < _source_x2) &&
            (y > _source_y1) && (y < _source_y2)) {
          // Add sample to metric
          *ptr2tmp =  round(_interpolator->EvaluateInside(x, y, 0));
          _metric->Add(*ptr2target, *ptr2tmp);
        } else {
          *ptr2tmp = -1;
        }
      }
      // Increment pointers to next voxel
      ptr2tmp++;
      ptr2target++;
      ptr += 3;
    }
  }

  // Evaluate similarity measure
  double similarity = _metric->Evaluate();

  // Add penalty for smoothness
  if (this->_Lambda1 > 0) {
    similarity += this->_Lambda1*this->SmoothnessPenalty();
  }
  // Add penalty for volume preservation
  if (this->_Lambda2 > 0) {
    similarity += this->_Lambda2*this->VolumePreservationPenalty();
  }
  // Add penalty for topology preservation
  if (this->_Lambda3 > 0) {
    similarity += this->_Lambda3*this->TopologyPreservationPenalty();
  }

  // Return similarity measure + penalty terms
  return similarity;
}

double irtkImageFreeFormRegistration2D::EvaluateDerivative(int index, double step)
{
  float *ptr;
  double bi, bj, bk, dx, dy, p[3];
  int i, j, i1, i2, j1, j2, k1, k2, dim;
  irtkPoint p1, p2;
  irtkGreyPixel *ptr2target, *ptr2tmp;

  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration2D::EvaluateDerivative(int, double)");

  // Initialize metrics for forward and backward derivative steps
  _tmpMetricA->ResetAndCopy(_metric);
  _tmpMetricB->ResetAndCopy(_metric);

  // Calculate bounding box of control point in world coordinates
  _affd->BoundingBoxCP(index, p1, p2);
  _target->WorldToImage(p1);
  _target->WorldToImage(p2);

  // Calculate bounding box of control point in image coordinates
  _affd->BoundingBoxImage(_target, index, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

  // Calculate incremental changes in lattice coordinates when looping
  // over target
  dx = (FFDLOOKUPTABLESIZE-1)/(p2._x-p1._x);
  dy = (FFDLOOKUPTABLESIZE-1)/(p2._y-p1._y);

  // Calculate whether this DOF corresponds to x, y or z-displacement
  dim = int(index / (_affd->GetX()*_affd->GetY()*_affd->GetZ()));

  // Loop over all voxels in the target (reference) volume
  bk = step * _localLookupTable[FFDLOOKUPTABLESIZE/2];
  for (j = j1; j <= j2; j++) {
    ptr2target = _target->GetPointerToVoxels(i1, j, 0);
    ptr        = &(_affdLookupTable[3*_target->VoxelToIndex(i1, j, 0)]);
    bj = bk * _localLookupTable[round(dy*(j-p1._y))];
    ptr2tmp  = _tmpImage->GetPointerToVoxels(i1, j, 0);
    for (i = i1; i <= i2; i++) {

      // Check whether reference point is valid
      if (*ptr2target >= 0) {
        bi = bj * _localLookupTable[round(dx*(i-p1._x))];

        // Delete old samples from both metrics
        if (*ptr2tmp != -1) {
          _tmpMetricA->Delete(*ptr2target, *ptr2tmp);
          _tmpMetricB->Delete(*ptr2target, *ptr2tmp);
        }

        p[0] = ptr[0];
        p[1] = ptr[1];
        p[2] = ptr[2];
        p[dim] += bi;

        // Convert transformed point to image coordinates
        _source->WorldToImage(p[0], p[1], p[2]);

        // Check whether transformed point is inside volume
        if ((p[0] > _source_x1) && (p[0] < _source_x2) &&
            (p[1] > _source_y1) && (p[1] < _source_y2)) {

          // Add sample to metric
          _tmpMetricA->Add(*ptr2target, round(_interpolator->EvaluateInside(p[0], p[1], 0)));
        }

        p[0] = ptr[0];
        p[1] = ptr[1];
        p[2] = ptr[2];
        p[dim] -= bi;

        // Convert transformed point to image coordinates
        _source->WorldToImage(p[0], p[1], p[2]);

        // Check whether transformed point is inside volume
        if ((p[0] > _source_x1) && (p[0] < _source_x2) &&
            (p[1] > _source_y1) && (p[1] < _source_y2)) {

          // Add sample to metric
          _tmpMetricB->Add(*ptr2target, round(_interpolator->EvaluateInside(p[0], p[1], 0)));
        }
      }

      // Increment pointers to next voxel
      ptr2target++;
      ptr2tmp++;
      ptr += 3;
    }
  }

  // Save value of DOF for which we calculate the derivative
  double dof = _affd->Get(index);

  // Evaluate similarity measure
  double similarityA = _tmpMetricA->Evaluate();

  // Add penalties
  _affd->Put(index, dof + step);

  // Smoothness
  if (this->_Lambda1 > 0) {
    similarityA += this->_Lambda1*this->SmoothnessPenalty(index);
  }
  // Volume preservation
  if (this->_Lambda2 > 0) {
    similarityA += this->_Lambda2*this->VolumePreservationPenalty(index);
  }
  // Topology preservation
  if (this->_Lambda3 > 0) {
    similarityA += this->_Lambda3*this->TopologyPreservationPenalty(index);
  }

  // Evaluate similarity measure
  double similarityB = _tmpMetricB->Evaluate();

  // Add penalties
  _affd->Put(index, dof - step);

  // Smoothness
  if (this->_Lambda1 > 0) {
    similarityB += this->_Lambda1*this->SmoothnessPenalty(index);
  }
  // Volume preservation
  if (this->_Lambda2 > 0) {
    similarityB += this->_Lambda2*this->VolumePreservationPenalty(index);
  }
  // Topology preservation
  if (this->_Lambda3 > 0) {
    similarityB += this->_Lambda3*this->TopologyPreservationPenalty(index);
  }

  // Restore value of DOF for which we calculate the derivative
  _affd->Put(index, dof);

  return similarityA - similarityB;
}
