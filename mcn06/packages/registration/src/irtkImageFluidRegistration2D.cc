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

// Used as temporary memory for transformed intensities
extern irtkGreyImage *_tmpSource;

// Used as temporary memory for metric
extern irtkSimilarityMetric *_tmpMetricA, *_tmpMetricB;

// Used as lookup table for transformed coordinates including level n. This
// lookup table needs to be calculated each time a control point has been
// modified.
extern float *_affdLookupTable;

// Used as lookup table for the contribution of each control point. This
// lookup table needs to be calculated only once.
extern float *_localLookupTable;

void irtkImageFluidRegistration2D::GuessParameter()
{
  cerr << "irtkImageFluidRegistration2D::GuessParameter: Not yet implemented" << endl;
  exit(1);
}

double irtkImageFluidRegistration2D::Evaluate()
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
  this->Debug("irtkImageFluidRegistration::Evaluate");

  // Initialize histogram
  _metric->Reset();

  // Loop over all voxels in the target (reference) volume
  ptr2target = _target->GetPointerToVoxels();
  ptr2tmp    = _tmpImage->GetPointerToVoxels();
  for (j = 0; j < _target->GetY(); j++) {
    for (i = 0; i < _target->GetX(); i++) {
      // Check whether reference point is valid
      if (*ptr2target >= 0) {
        x = i;
        y = j;
        z = 0;
        _target->ImageToWorld(x, y, z);
        _affd->Transform(x, y, z);
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

  // Return similarity measure + penalty terms
  return similarity;
}

double irtkImageFluidRegistration2D::EvaluateDerivative(int index, double step)
{
  float *ptr;
  double bi, bj, bk, dx, dy, p[3];
  int i, j, i1, i2, j1, j2, k1, k2, dim;
  irtkPoint p1, p2;
  irtkGreyPixel *ptr2target, *ptr2tmp;

  // Print debugging information
  this->Debug("irtkImageFluidRegistration2D::EvaluateDerivative(int, double)");

  // Initialize metrics for forward and backward derivative steps
  _tmpMetricA->Reset(_metric);
  _tmpMetricB->Reset(_metric);

  // Calculate bounding box of control point in world coordinates
  _affd->BoundingBox(index, p1, p2);
  _target->WorldToImage(p1);
  _target->WorldToImage(p2);

  // Calculate bounding box of control point in image coordinates
  _affd->BoundingBox(_target, index, i1, j1, k1, i2, j2, k2, this->_SpeedupFactor);

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

  // Evaluate similarity measure
  double similarityA = _tmpMetricA->Evaluate();

  // Evaluate similarity measure
  double similarityB = _tmpMetricB->Evaluate();

  return similarityA - similarityB;
}
