/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifdef HAS_TBB

extern concurrent_queue<irtkSimilarityMetric *> sim_queue;

class irtkMultiThreadedImageRigidRegistrationEvaluate2D
{

  /// Pointer to image transformation class
  irtkImageRigidRegistration *_filter;

  /// Pointer to metric
  irtkSimilarityMetric *_metric;

public:

  irtkMultiThreadedImageRigidRegistrationEvaluate2D(irtkImageRigidRegistration *filter) {
    // Initialize filter
    _filter = filter;

    // Initialize metric
    _metric = filter->_metric;
    _metric->Reset();
  }

  irtkMultiThreadedImageRigidRegistrationEvaluate2D(irtkMultiThreadedImageRigidRegistrationEvaluate2D &r, split) {

    // Initialize filter
    _filter = r._filter;

    // Copy similarity metric
    if (sim_queue.pop_if_present(_metric) == false) {
      _metric = irtkSimilarityMetric::New(_filter->_metric);
    }

    // Reset similarity metric
    _metric->Reset();
  }

  ~irtkMultiThreadedImageRigidRegistrationEvaluate2D() {
    if (_metric != _filter->_metric) sim_queue.push(_metric);
  }

  void join(irtkMultiThreadedImageRigidRegistrationEvaluate2D &rhs) {
    // Combine metrics
    _metric->Combine(rhs._metric);
  }

  void operator()(const blocked_range<int> &r) {
    int i, j;

    // Create iterator
    irtkHomogeneousTransformationIterator iterator((irtkHomogeneousTransformation *)_filter->_transformation);

    // Loop over all voxels in the target (reference) volume
    for (j = r.begin(); j != r.end(); j++) {

      // Initialize iterator
      iterator.Initialize(_filter->_target, _filter->_source, 0, j, 0);

      // Pointer to voxels in target image
      irtkGreyPixel *ptr2target = _filter->_target->GetPointerToVoxels(0, j, 0);

      // Loop over all voxels in the target (reference) volume
      for (j = 0; j < _filter->_target->GetY(); j++) {
        for (i = 0; i < _filter->_target->GetX(); i++) {
          // Check whether reference point is valid
          if (*ptr2target >= 0) {
            // Check whether transformed point is inside source volume
            if ((iterator._x > _filter->_source_x1) && (iterator._x < _filter->_source_x2) &&
                (iterator._y > _filter->_source_y1) && (iterator._y < _filter->_source_y2)) {
              // Add sample to metric
              _metric->Add(*ptr2target, round(_filter->_interpolator->EvaluateInside(iterator._x, iterator._y, 0)));
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
    }
  }
};

#endif
