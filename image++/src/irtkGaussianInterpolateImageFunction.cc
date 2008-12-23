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

#include <irtkGaussianInterpolateImageFunction.h>

template <class VoxelType>
irtkGaussianInterpolateImageFunction<VoxelType>::irtkGaussianInterpolateImageFunction(double Sigma)
{
  _Sigma = Sigma;
}

template <class VoxelType>
irtkGaussianInterpolateImageFunction<VoxelType>::~irtkGaussianInterpolateImageFunction(void)
{}

template <class VoxelType>
const char *irtkGaussianInterpolateImageFunction<VoxelType>::NameOfClass()
{
  return "irtkGaussianInterpolateImageFunction";
}

template <class VoxelType> void irtkGaussianInterpolateImageFunction<VoxelType>::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction<VoxelType>::Initialize();

  // Compute size of image
  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();
  this->_z = this->_input->GetZ();

  // Compute voxel size
  this->_input->GetPixelSize(&_xsize, &_ysize, &_zsize);

  // Compute filter extent
  _ExtentX = round(ceil(3*this->_Sigma/_xsize));
  _ExtentY = round(ceil(3*this->_Sigma/_ysize));
  _ExtentZ = round(ceil(3*this->_Sigma/_zsize));

  // Compute domain on which the Gaussian interpolation is defined
  this->_x1 = _ExtentX;
  this->_y1 = _ExtentY;
  this->_z1 = _ExtentZ;
  this->_x2 = this->_input->GetX() - _ExtentX;
  this->_y2 = this->_input->GetY() - _ExtentY;
  this->_z2 = this->_input->GetZ() - _ExtentZ;
}

template <class VoxelType> double irtkGaussianInterpolateImageFunction<VoxelType>::EvaluateInside(double x, double y, double z, double time)
{
  cerr << "irtkGaussianInterpolateImageFunction<VoxelType>::EvaluateInside: Not implemented" << endl;
  return 0;
}

template <class VoxelType> double irtkGaussianInterpolateImageFunction<VoxelType>::Evaluate(double x, double y, double z, double time)
{
  double val, sum;
  int x1, x2, y1, y2, z1, z2;

  // Initialize
  val = 0;
  sum = 0;
  x1 = x - _ExtentX;
  x2 = x + _ExtentX;
  y1 = y - _ExtentY;
  y2 = y + _ExtentY;
  z1 = z - _ExtentZ;
  z2 = z + _ExtentZ;

  // Create Gaussian kernel
  irtkScalarGaussian gaussian(this->_Sigma/_xsize, this->_Sigma/_ysize, this->_Sigma/_zsize, x, y, z);

  // Check whether boundary checking is necessary
  if ((x1 > 0) && (x2 < this->_input->GetX()) &&
      (y1 > 0) && (y2 < this->_input->GetY()) &&
      (z1 > 0) && (z2 < this->_input->GetZ())) {

    // If no, do fast convolution
    for (z = z1; z <= z2; z++) {
      for (y = y1; y <= y2; y++) {
        for (x = x1; x <= x2; x++) {
          double temp = gaussian.Evaluate(x, y, z);
          val += temp * this->_input->Get(x, y, z, time);
          sum += temp;
        }
      }
    }
  } else {
    // If yes, do slow convolution which handles boundaries
    for (z = z1; z <= z2; z++) {
      for (y = y1; y <= y2; y++) {
        for (x = x1; x <= x2; x++) {
          if ((x >= 0) && (x < this->_input->GetX()) &&
              (y >= 0) && (y < this->_input->GetY()) &&
              (z >= 0) && (z < this->_input->GetZ())) {
            double temp = gaussian.Evaluate(x, y, z);
            val += temp * this->_input->Get(x, y, z, time);
            sum += temp;
          }
        }
      }
    }
  }

  //  Normalize filter value by sum of filter elements
  if (sum > 0) {
    return val / sum;
  } else {
    return 0;
  }

}

template class irtkGaussianInterpolateImageFunction<irtkBytePixel>;
template class irtkGaussianInterpolateImageFunction<irtkGreyPixel>;
template class irtkGaussianInterpolateImageFunction<irtkRealPixel>;
