/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkMultiThreadedImageFreeFormRegistrationWithPadding.h 900 2013-05-29 16:26:00Z as12312 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2013-05-29 17:26:00 +0100 (Wed, 29 May 2013) $
  Version   : $Revision: 900 $
  Changes   : $Author: as12312 $

=========================================================================*/

#ifdef HAS_TBB

extern tbb::deprecated::concurrent_queue<irtkSimilarityMetric *> sim_queue;

class irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluate
{

  /// Pointer to image transformation class
  irtkImageFreeFormRegistrationWithPadding *_filter;

  /// Pointer to metric
  irtkSimilarityMetric *_metric;

public:

  irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluate(irtkImageFreeFormRegistrationWithPadding *filter) {
    // Initialize filter
    _filter = filter;

    // Initialize metric
    _metric = filter->_metric;
    _metric->Reset();
  }

  irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluate(irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluate &r, split) {

    // Initialize filter
    _filter = r._filter;

    // Copy similarity metric
    if (sim_queue.pop_if_present(_metric) == false) {
      _metric = irtkSimilarityMetric::New(_filter->_metric);
    }

    // Reset similarity metric
    _metric->Reset();
  }

  ~irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluate() {
    if (_metric != _filter->_metric) sim_queue.push(_metric);
  }

  void join(irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluate &rhs) {
    // Combine metrics
    _metric->Combine(rhs._metric);
  }

  void operator()(const blocked_range<int> &r) {
    // Image coordinates
    int i, j, k, t;
    // World coordinates
    double x, y, z;
    // Pointer to reference data
    irtkGreyPixel *ptr2target;
    irtkGreyPixel *ptr2tmp;
    float *ptr;

    // Loop over all voxels in the target (reference) volume
    for (t = 0; t < _filter->_target->GetT(); t++) {
      for (k = r.begin(); k != r.end(); k++) {
        ptr2target = _filter->_target->GetPointerToVoxels(0, 0, k, t);
        ptr2tmp    = _tmpImage->GetPointerToVoxels(0, 0, k, t);
        ptr        = &(_filter->_mffdLookupTable[3*k*_filter->_target->GetX()*_filter->_target->GetY()]);
        for (j = 0; j < _filter->_target->GetY(); j++) {
          for (i = 0; i < _filter->_target->GetX(); i++) {
            // Check whether reference point is valid
            if (*ptr2target >= 0) {
              x = i;
              y = j;
              z = k;
              _filter->_target->ImageToWorld(x, y, z);
              _filter->_affd->LocalDisplacement(x, y, z);
              x += ptr[0];
              y += ptr[1];
              z += ptr[2];
              _filter->_source->WorldToImage(x, y, z);
              // Check whether transformed point is inside volume
              if ((x > _filter->_source_x1) && (x < _filter->_source_x2) &&
                  (y > _filter->_source_y1) && (y < _filter->_source_y2) &&
                  (z > _filter->_source_z1) && (z < _filter->_source_z2)) {
                // Add sample to metric
                double value = _filter->_interpolator->EvaluateInside(x, y, z, t);
	        if (value>=0)
	        {
                  *ptr2tmp =  round(value);
                  _metric->Add(*ptr2target, *ptr2tmp);
  	        }
	        else 
		  *ptr2tmp = -1;
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
      }
    }
  }
};

class irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluateGradient
{

  /// Pointer to image transformation class
  irtkImageFreeFormRegistrationWithPadding *_filter;

  float _step, *_dx;

public:

  irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluateGradient(irtkImageFreeFormRegistrationWithPadding *filter, float *dx, float step) {
    _dx     = dx;
    _step   = step;
    _filter = filter;
  }

  void operator()(const blocked_range<int> &r) const {
    int i;

    // Loop over all voxels in the target (reference) volume
    for (i = r.begin(); i != r.end(); i++) {
      if (_filter->_affd->irtkTransformation::GetStatus(i) == _Active) {
        _dx[i] = _filter->EvaluateDerivative(i, _step);
      } else {
        _dx[i] = 0;
      }
    }
  }
};

#endif
