/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkHessianImageFilter.h>

template <class VoxelType> irtkHessianImageFilter<VoxelType>::irtkHessianImageFilter(int type)
{
  _type = type;

  _Padding = MIN_GREY;
}

template <class VoxelType> bool irtkHessianImageFilter<VoxelType>::RequiresBuffering(void)
{
  return true;
}

template <class VoxelType> const char *irtkHessianImageFilter<VoxelType>::NameOfClass()
{
  return "irtkHessianImageFilter";
}

template <class VoxelType> void irtkHessianImageFilter<VoxelType>::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageToImage::Initialize");

  // Check inputs and outputs
  if (this->_input == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input" << endl;
    exit(1);
  }

  if (this->_output == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no output" << endl;
    exit(1);
  }

  if (this->_input->IsEmpty() == true) {
    cerr << this->NameOfClass() << "::Run: Input is empty" << endl;
    exit(1);
  }

  if (this->_input->GetT() > 1) {
    cerr << this->NameOfClass() << "::Run: Only implemented for images with t = 1" << endl;
    exit(1);
  }

  // Check whether filter requires buffering
  if (this->RequiresBuffering()) {
    this->Debug("irtkHessianImageFilter::Initialize: Filter requires buffering");

    // Check whether filter has external buffer
    if (this->_input == this->_output) {
      this->Debug("irtkHessianImageFilter::Initialize: Filter has internal buffer");
      this->_tmp    = this->_output;
      this->_output = new irtkGenericImage<VoxelType>;
    } else {
      this->Debug("irtkHessianImageFilter::Initialize: Filter has external buffer");
      this->_tmp    = NULL;
    }
  } else {
    this->Debug("irtkHessianImageFilter::Initialize: Filter requires no buffering");
  }

  // Make sure that output has the correct dimensions
  if (_type == HESSIAN_VECTOR) {
    irtkImageAttributes attr = this->_input->GetImageAttributes();
    attr._t = 6;
    this->_output->Initialize(attr);
  } else {
    if (this->_input != this->_output) this->_output->Initialize(this->_input->GetImageAttributes());
  }
}

template <class VoxelType> void irtkHessianImageFilter<VoxelType>::Run()
{
  double dxx, dxy, dxz, dyy, dyz, dzz;
  int x, y, z, x1, y1, z1, x2, y2, z2;

  // Do the initial set up
  this->Initialize();

  for (z = 0; z < this->_input->GetZ(); ++z) {
    z1 = z - 1;
    if (z1 < 0) z1 = 0;
    z2 = z + 1;
    if (z2 > this->_input->GetZ()-1) z2 = this->_input->GetZ()-1;

    for (y = 0; y < this->_input->GetY(); ++y) {
      y1 = y - 1;
      if (y1 < 0) y1 = 0;
      y2 = y + 1;
      if (y2 > this->_input->GetY()-1) y2 = this->_input->GetY()-1;

      for (x = 0; x < this->_input->GetX(); ++x) {
        x1 = x - 1;
        if (x1 < 0) x1 = 0;
        x2 = x + 1;
        if (x2 > this->_input->GetX()-1) x2 = this->_input->GetX()-1;

        // Compute derivatives
        if (x1 != x2 &&
            this->_input->Get(x2, y, z) > _Padding &&
            this->_input->Get(x1, y, z) > _Padding) {
          dxx = (this->_input->Get(x2, y, z) - 2.0 * this->_input->Get(x, y, z) + this->_input->Get(x1, y, z)) / (this->_input->GetXSize() * this->_input->GetXSize());
        } else {
          dxx = 0;
        }
        if (x1 != x2 &&
        	y1 != y2 &&
			this->_input->Get(x2, y, z) > _Padding &&
			this->_input->Get(x1, y, z) > _Padding &&
			this->_input->Get(x, y2, z) > _Padding &&
			this->_input->Get(x, y1, z) > _Padding) {
		  dxy = (this->_input->Get(x2, y2, z) - this->_input->Get(x2, y1, z) - this->_input->Get(x1, y2, z) + this->_input->Get(x1, y1, z)) / ((x2 - x1) * (y2 - y1) * this->_input->GetXSize() * this->_input->GetYSize());
		} else {
		  dxy = 0;
		}
        if (x1 != x2 &&
			z1 != z2 &&
			this->_input->Get(x2, y, z) > _Padding &&
			this->_input->Get(x1, y, z) > _Padding &&
			this->_input->Get(x, y, z2) > _Padding &&
			this->_input->Get(x, y, z1) > _Padding) {
		  dxz = (this->_input->Get(x2, y, z2) - this->_input->Get(x2, y, z1) - this->_input->Get(x1, y, z2) + this->_input->Get(x1, y, z1)) / ((x2 - x1) * (z2 - z1) * this->_input->GetXSize() * this->_input->GetZSize());
		} else {
		  dxz = 0;
		}

        if (y1 != y2 &&
            this->_input->Get(x, y2, z) > _Padding &&
            this->_input->Get(x, y1, z) > _Padding) {
          dyy = (this->_input->Get(x, y2, z) - 2.0 * this->_input->Get(x, y, z) + this->_input->Get(x, y1, z)) / (this->_input->GetYSize() * this->_input->GetYSize());
        } else {
          dyy = 0;
        }
        if (y1 != y2 &&
			z1 != z2 &&
			this->_input->Get(x, y2, z) > _Padding &&
			this->_input->Get(x, y1, z) > _Padding &&
			this->_input->Get(x, y, z2) > _Padding &&
			this->_input->Get(x, y, z1) > _Padding) {
		  dyz = (this->_input->Get(x, y2, z2) - this->_input->Get(x, y2, z1) - this->_input->Get(x, y1, z2) + this->_input->Get(x, y1, z1)) / ((y2 - y1) * (z2 - z1) * this->_input->GetYSize() * this->_input->GetZSize());
		} else {
		  dyz = 0;
		}

        if (z1 != z2 &&
            this->_input->Get(x, y, z2) > _Padding &&
            this->_input->Get(x, y, z1) > _Padding) {
          dzz = (this->_input->Get(x, y, z2) - 2.0 * this->_input->Get(x, y, z) + this->_input->Get(x, y, z1)) / (this->_input->GetZSize() * this->_input->GetZSize());
        } else {
          dzz = 0;
        }

        switch (_type) {
          case HESSIAN_XX:
            this->_output->PutAsDouble(x, y, z, 0, dxx);
            break;
          case HESSIAN_XY:
            this->_output->PutAsDouble(x, y, z, 0, dxy);
            break;
          case HESSIAN_XZ:
            this->_output->PutAsDouble(x, y, z, 0, dxz);
            break;
          case HESSIAN_YY:
			this->_output->PutAsDouble(x, y, z, 0, dyy);
			break;
          case HESSIAN_YZ:
			this->_output->PutAsDouble(x, y, z, 0, dyz);
			break;
          case HESSIAN_ZZ:
			this->_output->PutAsDouble(x, y, z, 0, dzz);
			break;
          case HESSIAN_VECTOR:
        	this->_output->PutAsDouble(x, y, z, 0, dxx);
        	this->_output->PutAsDouble(x, y, z, 1, dxy);
        	this->_output->PutAsDouble(x, y, z, 2, dxz);
        	this->_output->PutAsDouble(x, y, z, 3, dyy);
        	this->_output->PutAsDouble(x, y, z, 4, dyz);
			this->_output->PutAsDouble(x, y, z, 5, dzz);
			break;
          default:
            cerr << this->NameOfClass() << "::Run: Unknown gradient computation" << endl;
            exit(1);
        }
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();
}

template class irtkHessianImageFilter<unsigned char>;
template class irtkHessianImageFilter<short>;
template class irtkHessianImageFilter<float>;
template class irtkHessianImageFilter<double>;
