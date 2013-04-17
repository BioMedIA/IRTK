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

class irtkMultiThreadedImageRigidRegistrationWithPaddingEvaluate
{

  /// Pointer to image transformation class
  irtkImageRigidRegistrationWithPadding *_filter;

  /// Pointer to metric
  irtkSimilarityMetric *_metric;

public:

  irtkMultiThreadedImageRigidRegistrationWithPaddingEvaluate(irtkImageRigidRegistrationWithPadding *filter) {
    // Initialize filter
    _filter = filter;

    // Initialize metric
    _metric = filter->_metric;
    _metric->Reset();
  }

  irtkMultiThreadedImageRigidRegistrationWithPaddingEvaluate(irtkMultiThreadedImageRigidRegistrationWithPaddingEvaluate &r, split) {

    // Initialize filter
    _filter = r._filter;

    // Copy similarity metric
    if (sim_queue.pop_if_present(_metric) == false) {
      _metric = irtkSimilarityMetric::New(_filter->_metric);
    }

    // Reset similarity metric
    _metric->Reset();
  }

  ~irtkMultiThreadedImageRigidRegistrationWithPaddingEvaluate() {
    if (_metric != _filter->_metric) sim_queue.push(_metric);
  }

  void join(irtkMultiThreadedImageRigidRegistrationWithPaddingEvaluate &rhs) {
    // Combine metrics
    _metric->Combine(rhs._metric);
  }

  void operator()(const blocked_range<int> &r) {
    int i, j, k, t;

    // Create iterator
    irtkHomogeneousTransformationIterator iterator((irtkHomogeneousTransformation *)_filter->_transformation);

    // Loop over all voxels in the target (reference) volume

    for (t = 0; t < _filter->_target->GetT(); t++) {
      for (k = r.begin(); k != r.end(); k++) {

        // Initialize iterator
        iterator.Initialize(_filter->_target, _filter->_source, 0, 0, k);

        // Pointer to voxels in target image
        irtkGreyPixel *ptr2target = _filter->_target->GetPointerToVoxels(0, 0, k, t);

        for (j = 0; j < _filter->_target->GetY(); j++) {
          for (i = 0; i < _filter->_target->GetX(); i++) {
            // Check whether reference point is valid
            if (*ptr2target >= 0) {
              // Check whether transformed point is inside source volume
              if ((iterator._x > _filter->_source_x1) && (iterator._x < _filter->_source_x2) &&
                  (iterator._y > _filter->_source_y1) && (iterator._y < _filter->_source_y2) &&
                  (iterator._z > _filter->_source_z1) && (iterator._z < _filter->_source_z2)) {
                // Add sample to metric
                //_metric->Add(*ptr2target, round(_filter->_interpolator->EvaluateInside(iterator._x, iterator._y, iterator._z, t)));
                double value = _filter->_interpolator->EvaluateInside(iterator._x, iterator._y, iterator._z, t);
	            if (value >= 0)
                  _metric->Add(*ptr2target, round(value));
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
  }
};

#endif
