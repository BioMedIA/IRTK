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
#include <irtkNoise.h>

template <class VoxelType> irtkRicianNoiseWithPadding<VoxelType>::irtkRicianNoiseWithPadding() : irtkRicianNoise<VoxelType>()
{
  _PaddingValue = std::numeric_limits<VoxelType>::min();
}

template <class VoxelType> irtkRicianNoiseWithPadding<VoxelType>::irtkRicianNoiseWithPadding(double Amplitude, VoxelType PaddingValue) : irtkRicianNoise<VoxelType>(Amplitude)
{
  _PaddingValue = PaddingValue;
}

template <class VoxelType> const char *irtkRicianNoiseWithPadding<VoxelType>::NameOfClass()
{
  return "irtkRicianNoiseWithPadding";
}

template <class VoxelType> double irtkRicianNoiseWithPadding<VoxelType>::Run(int x, int y, int z, int t)
{
  if (this->_input->Get(x, y, z, t) > this->_PaddingValue) {
    return this->irtkRicianNoise<VoxelType>::Run(x, y, z, t);
  } else {
    return this->_PaddingValue;
  }
}

template class irtkRicianNoiseWithPadding<irtkBytePixel>;
template class irtkRicianNoiseWithPadding<irtkGreyPixel>;
template class irtkRicianNoiseWithPadding<irtkRealPixel>;
