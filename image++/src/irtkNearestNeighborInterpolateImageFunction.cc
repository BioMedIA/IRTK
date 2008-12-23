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

#include <irtkImageFunction.h>

template <class VoxelType> irtkNearestNeighborInterpolateImageFunction<VoxelType>::irtkNearestNeighborInterpolateImageFunction()
{
}

template <class VoxelType> irtkNearestNeighborInterpolateImageFunction<VoxelType>::~irtkNearestNeighborInterpolateImageFunction(void)
{
}

template <class VoxelType> const char *irtkNearestNeighborInterpolateImageFunction<VoxelType>::NameOfClass()
{
  return "irtkNearestNeighborInterpolateImageFunction";
}

template <class VoxelType> void irtkNearestNeighborInterpolateImageFunction<VoxelType>::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction<VoxelType>::Initialize();

  // Compute image domain
  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();
  this->_z = this->_input->GetZ();

  // Compute domain on which the linear interpolation is defined
  this->_x1 = -0.5;
  this->_y1 = -0.5;
  this->_z1 = -0.5;
  this->_x2 = this->_input->GetX()-0.5;
  this->_y2 = this->_input->GetY()-0.5;
  this->_z2 = this->_input->GetZ()-0.5;
}

template <class VoxelType> double irtkNearestNeighborInterpolateImageFunction<VoxelType>::EvaluateInside(double x, double y, double z, double t)
{
  int i, j, k, l;

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  return this->_input->operator()(i, j, k, l);
}

template <class VoxelType> double irtkNearestNeighborInterpolateImageFunction<VoxelType>::Evaluate(double x, double y, double z, double t)
{
  int i, j, k, l;

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  if ((i < 0) || (i >= this->_x) || (j < 0) || (j >= this->_y) || (k < 0) || (k >= this->_z)) {
    return this->_DefaultValue;
  } else {
    return this->_input->operator()(i, j, k, l);
  }
}

template class irtkNearestNeighborInterpolateImageFunction<irtkBytePixel>;
template class irtkNearestNeighborInterpolateImageFunction<irtkGreyPixel>;
template class irtkNearestNeighborInterpolateImageFunction<irtkRealPixel>;
