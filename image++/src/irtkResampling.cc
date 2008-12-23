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

#include <irtkResampling.h>

#ifdef HAS_TBB

template <class VoxelType> class irtkMultiThreadedResampling
{

  /// Time frame to transform
  int _t;

  /// Pointer to image transformation class
  irtkResampling<VoxelType> *_filter;

public:

  irtkMultiThreadedResampling(irtkResampling<VoxelType> *filter, int t) {
    _t = t;
    _filter = filter;
  }

  void operator()(const blocked_range<int> &r) const {
    int i, j, k;
    double x, y, z;

    for (k = r.begin(); k != r.end(); k++) {
      for (j = 0; j < _filter->_output->GetY(); j++) {
        for (i = 0; i < _filter->_output->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          _filter->_output->ImageToWorld(x, y, z);
          _filter->_input ->WorldToImage(x, y, z);
          _filter->_output->PutAsDouble(i, j, k, _t, _filter->_Interpolator->Evaluate(x, y, z, _t));
        }
      }
    }
  }
};

#endif

template <class VoxelType> irtkResampling<VoxelType>::irtkResampling(double new_xsize, double new_ysize, double new_zsize)
{
  _XSize = new_xsize;
  _YSize = new_ysize;
  _ZSize = new_zsize;
  _Interpolator = NULL;
}

template <class VoxelType> Bool irtkResampling<VoxelType>::RequiresBuffering(void)
{
  return True;
}

template <class VoxelType> const char *irtkResampling<VoxelType>::NameOfClass()
{
  return "irtkResampling";
}

template <class VoxelType> void irtkResampling<VoxelType>::Initialize()
{
  int new_x, new_y, new_z;
  double xaxis[3], yaxis[3], zaxis[3];
  double new_xsize, new_ysize, new_zsize;
  double old_xsize, old_ysize, old_zsize;

  // Do the initial set up
  this->irtkImageToImage<VoxelType>::Initialize();

  // Set up interpolator
  if (_Interpolator == NULL) {
    cerr << "irtkResampling::Initialize: No interpolator found!" << endl;
    exit(1);
  }
  this->_Interpolator->SetInput(this->_input);
  this->_Interpolator->Initialize();

  // Determine the old dimensions of the image
  this->_input->GetPixelSize(&old_xsize, &old_ysize, &old_zsize);

  // Determine the new dimensions of the image
  new_x = int(this->_input->GetX() * old_xsize / _XSize);
  new_y = int(this->_input->GetY() * old_ysize / _YSize);
  new_z = int(this->_input->GetZ() * old_zsize / _ZSize);

  // Determine the new voxel dimensions
  if (new_x < 1) {
    new_x     =  1;
    new_xsize =  old_xsize;
  } else {
    new_xsize = this->_XSize;
  }
  if (new_y < 1) {
    new_y     =  1;
    new_ysize =  old_ysize;
  } else {
    new_ysize = this->_YSize;
  }
  if (new_z < 1) {
    new_z     =  1;
    new_zsize =  old_zsize;
  } else {
    new_zsize = this->_ZSize;
  }

  // Allocate new image
  *this->_output = irtkGenericImage<VoxelType>(new_x, new_y, new_z, this->_input->GetT());

  // Set new voxel size
  this->_output->PutPixelSize(new_xsize, new_ysize, new_zsize);

  // Set new orientation
  this->_input ->GetOrientation(xaxis, yaxis, zaxis);
  this->_output->PutOrientation(xaxis, yaxis, zaxis);

  // Set new origin
  this->_output->PutOrigin(this->_input->GetOrigin());
}

template <class VoxelType> void irtkResampling<VoxelType>::Run()
{
#ifdef HAS_TBB
  int l;
#else
  int i, j, k, l;
  double x, y, z;
#endif

  // Do the initial set up
  this->Initialize();

#ifdef HAS_TBB
  task_scheduler_init init(tbb_no_threads);

  tick_count t_start = tick_count::now();
#endif

  for (l = 0; l < this->_output->GetT(); l++) {

#ifdef HAS_TBB
    parallel_for(blocked_range<int>(0, this->_output->GetZ(), 1), irtkMultiThreadedResampling<VoxelType>(this, l));
#else

    for (k = 0; k < this->_output->GetZ(); k++) {
      for (j = 0; j < this->_output->GetY(); j++) {
        for (i = 0; i < this->_output->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          this->_output->ImageToWorld(x, y, z);
          this->_input ->WorldToImage(x, y, z);
          this->_output->PutAsDouble(i, j, k, l, this->_Interpolator->Evaluate(x, y, z, l));
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

template class irtkResampling<irtkBytePixel>;
template class irtkResampling<irtkGreyPixel>;
template class irtkResampling<irtkRealPixel>;
