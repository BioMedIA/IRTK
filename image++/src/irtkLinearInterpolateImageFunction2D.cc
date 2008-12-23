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

template <class VoxelType> irtkLinearInterpolateImageFunction2D<VoxelType>::irtkLinearInterpolateImageFunction2D()
{
}

template <class VoxelType> irtkLinearInterpolateImageFunction2D<VoxelType>::~irtkLinearInterpolateImageFunction2D(void)
{
}

template <class VoxelType> const char *irtkLinearInterpolateImageFunction2D<VoxelType>::NameOfClass()
{
  return "irtkLinearInterpolateImageFunction2D";
}

template <class VoxelType> void irtkLinearInterpolateImageFunction2D<VoxelType>::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction<VoxelType>::Initialize();

  // Check if image is 2D
  if (this->_input->GetZ() != 1) {
    cerr << "irtkLinearInterpolateImageFunction2D::Initialize(): ";
    cerr << "Input image is not 2D" << endl;
    exit(1);
  }

  // Compute size of image
  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();

  // Compute domain on which the linear interpolation is defined
  this->_x1 = 0;
  this->_y1 = 0;
  this->_x2 = this->_input->GetX() - 1;
  this->_y2 = this->_input->GetY() - 1;

  // Calculate offsets for fast pixel access
  this->_offset1 = 0;
  this->_offset2 = 1;
  this->_offset3 = this->_input->GetX();
  this->_offset4 = this->_input->GetX()+1;
}

template <class VoxelType> double irtkLinearInterpolateImageFunction2D<VoxelType>::EvaluateInside(double x, double y, double z, double time)
{
  int i, j;
  VoxelType *ptr;
  double t1, t2, u1, u2;

  // Calculated integer coordinates
  i  = int(x);
  j  = int(y);

  // Calculated fractional coordinates
  t1 = x - i;
  u1 = y - j;
  t2 = 1 - t1;
  u2 = 1 - u1;

  // Get pointer to data
  ptr = this->_input->GetPointerToVoxels(i, j, round(z), round(time));

  // Linear interpolation
  return (t1 * (u2 * ptr[this->_offset2] + u1 * ptr[this->_offset4]) +
          t2 * (u2 * ptr[this->_offset1] + u1 * ptr[this->_offset3]));
}

template <class VoxelType> double irtkLinearInterpolateImageFunction2D<VoxelType>::Evaluate(double x, double y, double z, double time)
{
  double val;
  int i, j, k, l, m, t;

  i = (int)floor(x);
  j = (int)floor(y);
  k = round(z);
  t = round(time);

  val = 0;
  for (l = i; l <= i+1; l++) {
    if ((l >= 0) && (l < this->_x)) {
      for (m = j; m <= j+1; m++) {
        if ((m >= 0) && (m < this->_y)) {
          val += (1 - fabs(l - x))*(1 - fabs(m - y))*this->_input->Get(l, m, k, t);
        }
      }
    }
  }
  return val;
}

template class irtkLinearInterpolateImageFunction2D<irtkBytePixel>;
template class irtkLinearInterpolateImageFunction2D<irtkGreyPixel>;
template class irtkLinearInterpolateImageFunction2D<irtkRealPixel>;
