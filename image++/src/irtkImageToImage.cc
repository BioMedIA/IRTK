/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkImageToImage.h>

#ifdef HAS_TBB

#include <tbb/tick_count.h>

template <class VoxelType> class irtkMultiThreadedImageToImage
{

  /// Time frame to transform
  int _t;

  /// Pointer to image transformation class
  irtkImageToImage<VoxelType> *_filter;

public:

  irtkMultiThreadedImageToImage(irtkImageToImage<VoxelType> *filter, int t) {
    _t = t;
    _filter = filter;
  }

  void operator()(const blocked_range<int> &r) const {
    int i, j, k;

    for (k = r.begin(); k != r.end(); k++) {
      for (j = 0; j < _filter->_input->GetY(); j++) {
        for (i = 0; i < _filter->_input->GetX(); i++) {
          _filter->_output->PutAsDouble(i, j, k, _t, _filter->Run(i, j, k, _t));
        }
      }
    }
  }
};

#endif

template <class VoxelType> irtkImageToImage<VoxelType>::irtkImageToImage()
{
  // Set in- and outputs
  _input  = NULL;
  _output = NULL;

  // Default parameters
  _DebugFlag = false;
}

template <class VoxelType> irtkImageToImage<VoxelType>::~irtkImageToImage()
{
  // Set in- and outputs
  _input  = NULL;
  _output = NULL;
}

template <class VoxelType> void irtkImageToImage<VoxelType>::SetInput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << "irtkImageToImage::SetInput: Input is not an image\n";
    exit(1);
  }
}

template <class VoxelType> void irtkImageToImage<VoxelType>::SetOutput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _output = image;
  } else {
    cerr << "irtkImageToImage::SetOutput: Output is not an image\n";
    exit(1);
  }
}

template <class VoxelType> void irtkImageToImage<VoxelType>::Debug(const char *message)
{
  if (_DebugFlag == true) cout << message << endl;
}

template <class VoxelType> void irtkImageToImage<VoxelType>::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageToImage::Initialize");

  // Check inputs and outputs
  if (_input == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input" << endl;
    exit(1);
  }

  if (_output == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no output" << endl;
    exit(1);
  }

  if (_input->IsEmpty() == true) {
    cerr << this->NameOfClass() << "::Run: Input is empty" << endl;
    exit(1);
  }

  // Check whether filter requires buffering
  if (this->RequiresBuffering()) {
    this->Debug("irtkImageToImage::Initialize: Filter requires buffering");

    // Check whether filter has external buffer
    if (_input == _output) {
      this->Debug("irtkImageToImage::Initialize: Filter has internal buffer");
      _tmp    = _output;
      _output = new irtkGenericImage<VoxelType>;
    } else {
      this->Debug("irtkImageToImage::Initialize: Filter has external buffer");
      _tmp    = NULL;
    }
  } else {
    this->Debug("irtkImageToImage::Initialize: Filter requires no buffering");
  }

  // Make sure that output has the correct dimensions
  if (_input != _output) _output->Initialize(_input->GetImageAttributes());
}

template <class VoxelType> void irtkImageToImage<VoxelType>::Finalize()
{
  // Print debugging information
  this->Debug("irtkImageToImage::Finalize");

  // Check whether filter requires buffering
  if (this->RequiresBuffering()) {
    this->Debug("irtkImageToImage::Finalize: Filter requires buffering");

    // Check whether filter has internal buffer
    if (_tmp != NULL) {
      this->Debug("irtkImageToImage::Finalize: Filter has internal buffer");
      // Copy buffer
      *_tmp = *_output;
      // Delete buffer
      delete _output;
      // Bend pointers back
      _output = _tmp;
    } else {
      this->Debug("irtkImageToImage::Finalize: Filter has external buffer");
    }
  } else {
    this->Debug("irtkImageToImage::Finalize: Filter requires no buffering");
  }
}

template <class VoxelType> double irtkImageToImage<VoxelType>::Run(int, int, int, int)
{
  cerr << "Filter " << this->NameOfClass() << " has no Run(int, int, int) ";
  cerr << "member function. Using irtkImageToImage::Run." << endl;
  return 0;
}

template <class VoxelType> void irtkImageToImage<VoxelType>::Run()
{
#ifdef HAS_TBB
  int t;
#else
  int x, y, z, t;
#endif

  // Do the initial set up
  this->Initialize();

#ifdef HAS_TBB
  task_scheduler_init init(tbb_no_threads);

  tick_count t_start = tick_count::now();
#endif

  // Calculate
  for (t = 0; t < _input->GetT(); t++) {

#ifdef HAS_TBB
    parallel_for(blocked_range<int>(0, this->_output->GetZ(), 1), irtkMultiThreadedImageToImage<VoxelType>(this, t));
#else

    for (z = 0; z < _input->GetZ(); z++) {
      for (y = 0; y < _input->GetY(); y++) {
        for (x = 0; x < _input->GetX(); x++) {
          _output->PutAsDouble(x, y, z, t, this->Run(x, y, z, t));
        }
      }
    }

#endif

  }

#ifdef HAS_TBB

  tick_count t_end = tick_count::now();
  if (tbb_debug) cout << this->NameOfClass() << " = " << (t_end - t_start).seconds() << " secs." << endl;
  init.terminate();

#endif

  // Do the final cleaning up
  this->Finalize();
}

template class irtkImageToImage<unsigned char>;
template class irtkImageToImage<short>;
template class irtkImageToImage<unsigned short>;
template class irtkImageToImage<float>;
template class irtkImageToImage<double>;
