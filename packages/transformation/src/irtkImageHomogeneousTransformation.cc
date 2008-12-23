/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkTransformation.h>

#include <irtkHomogeneousTransformationIterator.h>

#ifdef HAS_TBB

#include <tbb/tick_count.h>

template <class VoxelType> class irtkMultiThreadedImageHomogeneousTransformation
{

  /// Time frame to transform
  int _toutput, _tinput;

  /// Pointer to image transformation class
  irtkImageHomogeneousTransformation<VoxelType> *_imagetransformation;

public:

  irtkMultiThreadedImageHomogeneousTransformation(irtkImageHomogeneousTransformation<VoxelType> *imagetransformation, int toutput, int tinput) {
    _toutput = toutput;
    _tinput  = tinput;
    _imagetransformation = imagetransformation;
  }

  void operator()(const blocked_range<int> &r) const {
    int i, j, k;
    irtkPoint p;

    for (k = r.begin(); k != r.end(); k++) {
      // Create iterator
      irtkHomogeneousTransformationIterator
      iterator((irtkHomogeneousTransformation *)_imagetransformation->_transformation);
      iterator.Initialize(_imagetransformation->_output, _imagetransformation->_input, 0, 0, k);

      // Pointer to voxels in output image
      VoxelType *ptr2output = _imagetransformation->_output->GetPointerToVoxels(0, 0, k, _toutput);

      for (j = 0; j < _imagetransformation->_output->GetY(); j++) {
        for (i = 0; i < _imagetransformation->_output->GetX(); i++) {
          if (*ptr2output > _imagetransformation->_TargetPaddingValue) {
            // Check whether transformed point is inside input volume
            if ((iterator._x > -0.5) && (iterator._x < _imagetransformation->_input->GetX()-0.5) &&
                (iterator._y > -0.5) && (iterator._y < _imagetransformation->_input->GetY()-0.5) &&
                (iterator._z > -0.5) && (iterator._z < _imagetransformation->_input->GetZ()-0.5)) {
              *ptr2output = round(_imagetransformation->_interpolator->Evaluate(iterator._x, iterator._y, iterator._z, _tinput));
            } else {
              *ptr2output = _imagetransformation->_SourcePaddingValue;
            }
          } else {
            *ptr2output = _imagetransformation->_SourcePaddingValue;
          }
          iterator.NextX();
          ptr2output++;
        }
        iterator.NextY();
      }
      iterator.NextZ();
    }
  }
};

#endif

template <class VoxelType> irtkImageHomogeneousTransformation<VoxelType>::irtkImageHomogeneousTransformation() : irtkImageTransformation<VoxelType>()
{}

template <class VoxelType> irtkImageHomogeneousTransformation<VoxelType>::~irtkImageHomogeneousTransformation()
{}

template <class VoxelType> void irtkImageHomogeneousTransformation<VoxelType>::SetTransformation(irtkTransformation *transformation)
{
  if (strcmp(transformation->NameOfClass(), "irtkRigidTransformation") == 0) {
    this->_transformation = (irtkHomogeneousTransformation *)transformation;
    return;
  }
  if (strcmp(transformation->NameOfClass(), "irtkAffineTransformation") == 0) {
    this->_transformation = (irtkHomogeneousTransformation *)transformation;
    return;
  }
  if (strcmp(transformation->NameOfClass(), "irtkHomogeneousTransformation") == 0) {
    this->_transformation = (irtkHomogeneousTransformation *)transformation;
    return;
  }
  cerr << "irtkImageHomogeneousTransformation::SetTransform: Not a ";
  cerr << "homogeneous transformation: " << transformation->NameOfClass()
       << endl;
  exit(1);
}

template <class VoxelType> void irtkImageHomogeneousTransformation<VoxelType>::Run()
{
  int i, j, k, l;
  VoxelType *ptr2output;

  // Check inputs and outputs
  if (this->_input == NULL) {
    cerr << "irtkImageHomogeneousTransformation::Run: Filter has no input" << endl;
    exit(1);
  }

  if (this->_output == NULL) {
    cerr << "irtkImageHomogeneousTransformation::Run: Filter has no output" << endl;
    exit(1);
  }

  if (this->_transformation == NULL) {
    cerr << "irtkImageHomogeneousTransformation::Run: Filter has no transformation" << endl;
    exit(1);
  }

  if (this->_interpolator == NULL) {
    cerr << "irtkImageHomogeneousTransformation::Run: Filter has no interpolator" << endl;
    exit(1);
  }

  if (this->_input->IsEmpty() == True) {
    cerr << "irtkImageHomogeneousTransformation::Run: Input is empty" << endl;
    exit(1);
  }

  if (this->_input == this->_output) {
    cerr << "irtkImageHomogeneousTransformation::Run: Input equals output" << endl;
    exit(1);
  }

  // Setup interpolation
  this->_interpolator->SetInput(this->_input);
  this->_interpolator->Initialize();

  // Invert transformation
  if (this->_Invert == True) ((irtkHomogeneousTransformation *)this->_transformation)->Invert();

  // Create iterator
  irtkHomogeneousTransformationIterator
  iterator((irtkHomogeneousTransformation *)this->_transformation);

  // Pointer to voxels in output image
  ptr2output = this->_output->GetPointerToVoxels();

#ifdef HAS_TBB
  task_scheduler_init init(tbb_no_threads);

  tick_count t_start = tick_count::now();
#endif

  // Loop over all voxels in the output (reference) volume
  for (l = 0; l < this->_output->GetT(); l++) {
    int t = round(this->_input->TimeToImage(this->_output->ImageToTime(l)));
    if ((t >= 0) || (t < this->_input->GetT())) {

#ifdef HAS_TBB
      parallel_for(blocked_range<int>(0, this->_output->GetZ(), 1), irtkMultiThreadedImageHomogeneousTransformation<VoxelType>(this, l, t));
#else

      // Initialize iterator
      iterator.Initialize(this->_output, this->_input);
      for (k = 0; k < this->_output->GetZ(); k++) {
        for (j = 0; j < this->_output->GetY(); j++) {
          for (i = 0; i < this->_output->GetX(); i++) {
            if (*ptr2output > this->_TargetPaddingValue) {
              // Check whether transformed point is inside input volume
              if ((iterator._x > -0.5) && (iterator._x < this->_input->GetX()-0.5) &&
                  (iterator._y > -0.5) && (iterator._y < this->_input->GetY()-0.5) &&
                  (iterator._z > -0.5) && (iterator._z < this->_input->GetZ()-0.5)) {
                *ptr2output = round(this->_interpolator->Evaluate(iterator._x, iterator._y, iterator._z, t));
              } else {
                *ptr2output = this->_SourcePaddingValue;
              }
            } else {
              *ptr2output = this->_SourcePaddingValue;
            }
            iterator.NextX();
            ptr2output++;
          }
          iterator.NextY();
        }
        iterator.NextZ();
      }

#endif

    } else {
      for (k = 0; k < this->_output->GetZ(); k++) {
        for (j = 0; j < this->_output->GetY(); j++) {
          for (i = 0; i < this->_output->GetX(); i++) {
            *ptr2output = this->_SourcePaddingValue;
            ptr2output++;
          }
        }
      }
    }
  }

#ifdef HAS_TBB

  tick_count t_end = tick_count::now();
  if (tbb_debug) cout << "irtkImageHomogeneousTransformation = " << (t_end - t_start).seconds() << " secs." << endl;
  init.terminate();

#endif

  // Invert transformation
  if (this->_Invert == True) ((irtkHomogeneousTransformation *)this->_transformation)->Invert();
}

template class irtkImageHomogeneousTransformation<irtkBytePixel>;
template class irtkImageHomogeneousTransformation<irtkGreyPixel>;
template class irtkImageHomogeneousTransformation<irtkRealPixel>;
