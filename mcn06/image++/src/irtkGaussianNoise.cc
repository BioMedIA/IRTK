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

#include <nr.h>
#include <nrutil.h>

template <class VoxelType> irtkGaussianNoise<VoxelType>::irtkGaussianNoise() : irtkNoise<VoxelType>()
{
  _Mean   = 0;
  _Sigma  = 1;
  _MinVal = VoxelType(MIN_GREY);
  _MaxVal = VoxelType(MAX_GREY);

  long temp = -1 * this->_Init;
  (void) gasdev(&temp);
}

template <class VoxelType> irtkGaussianNoise<VoxelType>::irtkGaussianNoise(double Mean, double Sigma, VoxelType MinVal, VoxelType MaxVal) : irtkNoise<VoxelType>()
{
  this->_Mean   = Mean;
  this->_Sigma  = Sigma;
  this->_MinVal = MinVal;
  this->_MaxVal = MaxVal;

  long temp = -1 * this->_Init;
  (void) gasdev(&temp);
}


template <class VoxelType> const char *irtkGaussianNoise<VoxelType>::NameOfClass()
{
  return "irtkGaussianNoise";
}

template <class VoxelType> double irtkGaussianNoise<VoxelType>::Run(int x, int y, int z, int t)
{
  double tmp = this->_input->Get(x, y, z, t) + this->_Sigma * gasdev(&this->_Init) + this->_Mean;
  if (tmp < this->_MinVal) return this->_MinVal;
  if (tmp > this->_MaxVal) return this->_MaxVal;
  return tmp;
}

template class irtkGaussianNoise<irtkBytePixel>;
template class irtkGaussianNoise<irtkGreyPixel>;
template class irtkGaussianNoise<irtkRealPixel>;

