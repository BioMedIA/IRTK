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

#include <irtkNoise.h>

template <class VoxelType> irtkGaussianNoiseWithPadding<VoxelType>::irtkGaussianNoiseWithPadding() : irtkGaussianNoise<VoxelType>()
{
  _PaddingValue = std::numeric_limits<VoxelType>::min();
}

template <class VoxelType> irtkGaussianNoiseWithPadding<VoxelType>::irtkGaussianNoiseWithPadding(double Mean, double Sigma, VoxelType PaddingValue, VoxelType MinVal, VoxelType MaxVal) : irtkGaussianNoise<VoxelType>(Mean, Sigma, MinVal, MaxVal)
{
  _PaddingValue = PaddingValue;
}

template <class VoxelType> const char *irtkGaussianNoiseWithPadding<VoxelType>::NameOfClass()
{
  return "irtkGaussianNoiseWithPadding";
}

template <class VoxelType> double irtkGaussianNoiseWithPadding<VoxelType>::Run(int x, int y, int z, int t)
{
  if (this->_input->Get(x, y, z, t) > this->_PaddingValue) {
    return this->irtkGaussianNoise<VoxelType>::Run(x, y, z, t);
  } else {
    return this->_PaddingValue;
  }
}

template class irtkGaussianNoiseWithPadding<irtkBytePixel>;
template class irtkGaussianNoiseWithPadding<irtkGreyPixel>;
template class irtkGaussianNoiseWithPadding<irtkRealPixel>;
