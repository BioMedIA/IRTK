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

template <class VoxelType> class irtkMultiThreadedResamplingWithPadding
{

  /// Time frame to transform
  int _t;

  /// Pointer to image transformation class
  irtkResamplingWithPadding<VoxelType> *_filter;

public:

  irtkMultiThreadedResamplingWithPadding(irtkResamplingWithPadding<VoxelType> *filter, int t) {
    _t = t;
    _filter = filter;
  }

  void operator()(const blocked_range<int> &r) const {
    int i, j, k, l, u, v, w, pad;
    double val, sum;
    double w1, w2, w3, w4, w5, w6, w7, w8, dx, dy, dz, x, y, z;

    l = _t;
    for (k = r.begin(); k != r.end(); k++) {
      for (j = 0; j < _filter->_output->GetY(); j++) {
        for (i = 0; i < _filter->_output->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          _filter->_output->ImageToWorld(x, y, z);
          _filter->_input ->WorldToImage(x, y, z);

          // Calculate integer fraction of points
          u = (int)floor(x);
          v = (int)floor(y);
          w = (int)floor(z);

          // Calculate floating point fraction of points
          dx = x - u;
          dy = y - v;
          dz = z - w;

          // Calculate weights for trilinear interpolation
          w1 = (1 - dx) * (1 - dy) * (1 - dz);
          w2 = (1 - dx) * (1 - dy) * dz;
          w3 = (1 - dx) * dy * (1 - dz);
          w4 = (1 - dx) * dy * dz;
          w5 = dx * (1 - dy) * (1 - dz);
          w6 = dx * (1 - dy) * dz;
          w7 = dx * dy * (1 - dz);
          w8 = dx * dy * dz;

          // Calculate trilinear interpolation, ignoring padded values
          val = 0;
          pad = 8;
          sum = 0;
          if ((u >= 0) && (u < _filter->_input->GetX()) &&
              (v >= 0) && (v < _filter->_input->GetY()) &&
              (w >= 0) && (w < _filter->_input->GetZ())) {
            if (_filter->_input->Get(u, v, w, l) != _filter->_PaddingValue) {
              pad--;
              val += _filter->_input->Get(u, v, w, l) * w1;
              sum += w1;
            }
          } else {
            pad--;
          }
          if ((u >= 0)   && (u < _filter->_input->GetX()) &&
              (v >= 0)   && (v < _filter->_input->GetY()) &&
              (w+1 >= 0) && (w+1 < _filter->_input->GetZ())) {
            if (_filter->_input->Get(u, v, w+1, l) != _filter->_PaddingValue) {
              pad--;
              val += _filter->_input->Get(u, v, w+1, l) * w2;
              sum += w2;
            }
          } else {
            pad--;
          }
          if ((u >= 0)   && (u < _filter->_input->GetX()) &&
              (v+1 >= 0) && (v+1 < _filter->_input->GetY()) &&
              (w >= 0)   && (w < _filter->_input->GetZ())) {
            if (_filter->_input->Get(u, v+1, w, l) != _filter->_PaddingValue) {
              pad--;
              val += _filter->_input->Get(u, v+1, w, l) * w3;
              sum += w3;
            }
          } else {
            pad--;
          }
          if ((u >= 0)   && (u < _filter->_input->GetX()) &&
              (v+1 >= 0) && (v+1 < _filter->_input->GetY()) &&
              (w+1 >= 0) && (w+1 < _filter->_input->GetZ())) {
            if (_filter->_input->Get(u, v+1, w+1, l) != _filter->_PaddingValue) {
              pad--;
              val += _filter->_input->Get(u, v+1, w+1, l) * w4;
              sum += w4;
            }
          } else {
            pad--;
          }
          if ((u+1 >= 0) && (u+1 < _filter->_input->GetX()) &&
              (v >= 0)   && (v < _filter->_input->GetY()) &&
              (w >= 0)   && (w < _filter->_input->GetZ())) {
            if (_filter->_input->Get(u+1, v, w, l) != _filter->_PaddingValue) {
              pad--;
              val += _filter->_input->Get(u+1, v, w, l) * w5;
              sum += w5;
            }
          } else {
            pad--;
          }
          if ((u+1 >= 0) && (u+1 < _filter->_input->GetX()) &&
              (v >= 0)   && (v < _filter->_input->GetY()) &&
              (w+1 >= 0) && (w+1 < _filter->_input->GetZ())) {
            if (_filter->_input->Get(u+1, v, w+1, l) != _filter->_PaddingValue) {
              pad--;
              val += _filter->_input->Get(u+1, v, w+1, l) * w6;
              sum += w6;
            }
          } else {
            pad--;
          }
          if ((u+1 >= 0) && (u+1 < _filter->_input->GetX()) &&
              (v+1 >= 0) && (v+1 < _filter->_input->GetY()) &&
              (w >= 0)   && (w < _filter->_input->GetZ())) {
            if (_filter->_input->Get(u+1, v+1, w, l) != _filter->_PaddingValue) {
              pad--;
              val += _filter->_input->Get(u+1, v+1, w, l) * w7;
              sum += w7;
            }
          } else {
            pad--;
          }
          if ((u+1 >= 0) && (u+1 < _filter->_input->GetX()) &&
              (v+1 >= 0) && (v+1 < _filter->_input->GetY()) &&
              (w+1 >= 0) && (w+1 < _filter->_input->GetZ())) {
            if (_filter->_input->Get(u+1, v+1, w+1, l) != _filter->_PaddingValue) {
              pad--;
              val += _filter->_input->Get(u+1, v+1, w+1, l) * w8;
              sum += w8;
            }
          } else {
            pad--;
          }
          if (pad < 4) {
            if (sum > 0) {
              _filter->_output->PutAsDouble(i, j, k, l, val / sum);
            } else {
              // Special case: Point lies on edge of padded voxels.
              _filter->_output->Put(i, j, k, l, _filter->_PaddingValue);
              // Rather than:
              // cerr << _filter->NameOfClass();
              // cerr << "::Run: _filter should never happen" << endl;
              // exit(1);
            }
          } else {
            _filter->_output->Put(i, j, k, l, _filter->_PaddingValue);
          }
        }
      }
    }
  }
};

#endif

template <class VoxelType>
irtkResamplingWithPadding<VoxelType>::irtkResamplingWithPadding(double new_xsize, double new_ysize, double new_zsize, VoxelType PaddingValue) : irtkResampling<VoxelType>(new_xsize, new_ysize, new_zsize)
{
  _PaddingValue = PaddingValue;
}

template <class VoxelType>
bool irtkResamplingWithPadding<VoxelType>::RequiresBuffering(void)
{
  return true;
}

template <class VoxelType>
const char *irtkResamplingWithPadding<VoxelType>::NameOfClass()
{
  return "irtkResamplingWithPadding";
}

template <class VoxelType>
void irtkResamplingWithPadding<VoxelType>::Initialize()
{
  int new_x, new_y, new_z;
  double xaxis[3], yaxis[3], zaxis[3];
  double new_xsize, new_ysize, new_zsize;
  double old_xsize, old_ysize, old_zsize;

  // Do the initial set up
  this->irtkImageToImage<VoxelType>::Initialize();

  // Determine the old dimensions of the image
  this->_input->GetPixelSize(&old_xsize, &old_ysize, &old_zsize);

  // Determine the new dimensions of the image
  new_x = int(this->_input->GetX() * old_xsize / this->_XSize);
  new_y = int(this->_input->GetY() * old_ysize / this->_YSize);
  new_z = int(this->_input->GetZ() * old_zsize / this->_ZSize);

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

template <class VoxelType> void irtkResamplingWithPadding<VoxelType>::Run()
{
#ifdef HAS_TBB
  int l;
#else
  int i, j, k, l, u, v, w, pad;
  double val, sum;
  double w1, w2, w3, w4, w5, w6, w7, w8, dx, dy, dz, x, y, z;
#endif

  // Do the initial set up
  this->Initialize();

#ifdef HAS_TBB
  task_scheduler_init init(tbb_no_threads);

  tick_count t_start = tick_count::now();
#endif

  for (l = 0; l < this->_output->GetT(); l++) {

#ifdef HAS_TBB
    parallel_for(blocked_range<int>(0, this->_output->GetZ(), 1), irtkMultiThreadedResamplingWithPadding<VoxelType>(this, l));
#else

    for (k = 0; k < this->_output->GetZ(); k++) {
      for (j = 0; j < this->_output->GetY(); j++) {
        for (i = 0; i < this->_output->GetX(); i++) {
          x = i;
          y = j;
          z = k;
          this->_output->ImageToWorld(x, y, z);
          this->_input ->WorldToImage(x, y, z);

          // Calculate integer fraction of points
          u = (int)floor(x);
          v = (int)floor(y);
          w = (int)floor(z);

          // Calculate floating point fraction of points
          dx = x - u;
          dy = y - v;
          dz = z - w;

          // Calculate weights for trilinear interpolation
          w1 = (1 - dx) * (1 - dy) * (1 - dz);
          w2 = (1 - dx) * (1 - dy) * dz;
          w3 = (1 - dx) * dy * (1 - dz);
          w4 = (1 - dx) * dy * dz;
          w5 = dx * (1 - dy) * (1 - dz);
          w6 = dx * (1 - dy) * dz;
          w7 = dx * dy * (1 - dz);
          w8 = dx * dy * dz;

          // Calculate trilinear interpolation, ignoring padded values
          val = 0;
          pad = 8;
          sum = 0;
          if ((u >= 0) && (u < this->_input->GetX()) &&
              (v >= 0) && (v < this->_input->GetY()) &&
              (w >= 0) && (w < this->_input->GetZ())) {
            if (this->_input->Get(u, v, w, l) != this->_PaddingValue) {
              pad--;
              val += this->_input->Get(u, v, w, l) * w1;
              sum += w1;
            }
          } else {
            pad--;
          }
          if ((u >= 0)   && (u < this->_input->GetX()) &&
              (v >= 0)   && (v < this->_input->GetY()) &&
              (w+1 >= 0) && (w+1 < this->_input->GetZ())) {
            if (this->_input->Get(u, v, w+1, l) != this->_PaddingValue) {
              pad--;
              val += this->_input->Get(u, v, w+1, l) * w2;
              sum += w2;
            }
          } else {
            pad--;
          }
          if ((u >= 0)   && (u < this->_input->GetX()) &&
              (v+1 >= 0) && (v+1 < this->_input->GetY()) &&
              (w >= 0)   && (w < this->_input->GetZ())) {
            if (this->_input->Get(u, v+1, w, l) != this->_PaddingValue) {
              pad--;
              val += this->_input->Get(u, v+1, w, l) * w3;
              sum += w3;
            }
          } else {
            pad--;
          }
          if ((u >= 0)   && (u < this->_input->GetX()) &&
              (v+1 >= 0) && (v+1 < this->_input->GetY()) &&
              (w+1 >= 0) && (w+1 < this->_input->GetZ())) {
            if (this->_input->Get(u, v+1, w+1, l) != this->_PaddingValue) {
              pad--;
              val += this->_input->Get(u, v+1, w+1, l) * w4;
              sum += w4;
            }
          } else {
            pad--;
          }
          if ((u+1 >= 0) && (u+1 < this->_input->GetX()) &&
              (v >= 0)   && (v < this->_input->GetY()) &&
              (w >= 0)   && (w < this->_input->GetZ())) {
            if (this->_input->Get(u+1, v, w, l) != this->_PaddingValue) {
              pad--;
              val += this->_input->Get(u+1, v, w, l) * w5;
              sum += w5;
            }
          } else {
            pad--;
          }
          if ((u+1 >= 0) && (u+1 < this->_input->GetX()) &&
              (v >= 0)   && (v < this->_input->GetY()) &&
              (w+1 >= 0) && (w+1 < this->_input->GetZ())) {
            if (this->_input->Get(u+1, v, w+1, l) != this->_PaddingValue) {
              pad--;
              val += this->_input->Get(u+1, v, w+1, l) * w6;
              sum += w6;
            }
          } else {
            pad--;
          }
          if ((u+1 >= 0) && (u+1 < this->_input->GetX()) &&
              (v+1 >= 0) && (v+1 < this->_input->GetY()) &&
              (w >= 0)   && (w < this->_input->GetZ())) {
            if (this->_input->Get(u+1, v+1, w, l) != this->_PaddingValue) {
              pad--;
              val += this->_input->Get(u+1, v+1, w, l) * w7;
              sum += w7;
            }
          } else {
            pad--;
          }
          if ((u+1 >= 0) && (u+1 < this->_input->GetX()) &&
              (v+1 >= 0) && (v+1 < this->_input->GetY()) &&
              (w+1 >= 0) && (w+1 < this->_input->GetZ())) {
            if (this->_input->Get(u+1, v+1, w+1, l) != this->_PaddingValue) {
              pad--;
              val += this->_input->Get(u+1, v+1, w+1, l) * w8;
              sum += w8;
            }
          } else {
            pad--;
          }
          if (pad < 4) {
            if (sum > 0) {
              this->_output->PutAsDouble(i, j, k, l, val / sum);
            } else {
              // Special case: Point lies on edge of padded voxels.
              this->_output->Put(i, j, k, l, this->_PaddingValue);
              // Rather than:
              // cerr << this->NameOfClass();
              // cerr << "::Run: This should never happen" << endl;
              // exit(1);
            }
          } else {
            this->_output->Put(i, j, k, l, this->_PaddingValue);
          }
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

template class irtkResamplingWithPadding<char>;
template class irtkResamplingWithPadding<unsigned char>;
template class irtkResamplingWithPadding<short>;
template class irtkResamplingWithPadding<unsigned short>;
template class irtkResamplingWithPadding<int>;
template class irtkResamplingWithPadding<float>;
template class irtkResamplingWithPadding<double>;
