/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkImage.h>

#include <irtkGradientImage.h>

template <class VoxelType> irtkGradientImage<VoxelType>::irtkGradientImage()
{
  _Padding = MIN_GREY;
}

template <class VoxelType> bool irtkGradientImage<VoxelType>::RequiresBuffering(void)
{
  return true;
}

template <class VoxelType> const char *irtkGradientImage<VoxelType>::NameOfClass()
{
  return "irtkGradientImage";
}

template <class VoxelType> double irtkGradientImage<VoxelType>::Run(int x, int y, int z, int t)
{
  double dx, dy, dz;

  if ((x > 0) && (x < this->_input->GetX()-1)
      && this->_input->Get(x-1, y, z, t) > _Padding
      && this->_input->Get(x+1, y, z, t) > _Padding
     ) {
    dx = this->_input->Get(x-1, y, z, t) - this->_input->Get(x+1, y, z, t);
  } else {
    dx = 0;
  }

  if ((y > 0) && (y < this->_input->GetY()-1)
      && this->_input->Get(x, y-1, z, t) > _Padding
      && this->_input->Get(x, y+1, z, t) > _Padding
     ) {
    dy = this->_input->Get(x, y-1, z, t) - this->_input->Get(x, y+1, z, t);
  } else {
    dy = 0;
  }

  if ((z > 0) && (z < this->_input->GetZ()-1)
      && this->_input->Get(x, y, z-1, t) > _Padding
      && this->_input->Get(x, y, z+1, t) > _Padding
     ) {
    dz = this->_input->Get(x, y, z-1, t) - this->_input->Get(x, y, z+1, t);
  } else {
    dz = 0;
  }

  return sqrt(dx*dx + dy*dy + dz*dz);
}

template <class VoxelType> void irtkGradientImage<VoxelType>::Run()
{
  int x, y, z, t;

  // Do the initial set up
  this->Initialize();

  for (t = 0; t < this->_input->GetT(); ++t) {
    for (z = 0; z < this->_input->GetZ(); ++z) {
      for (y = 0; y < this->_input->GetY(); ++y) {
        for (x = 0; x < this->_input->GetX(); ++x) {
          this->_output->PutAsDouble(x, y, z, t, this->Run(x, y, z, t));
        }
      }
    }
  }
  // Do the final cleaning up
  this->Finalize();
}

template class irtkGradientImage<unsigned char>;
template class irtkGradientImage<short>;
template class irtkGradientImage<unsigned short>;
template class irtkGradientImage<float>;
template class irtkGradientImage<double>;
