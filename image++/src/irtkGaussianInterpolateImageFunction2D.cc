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

#include <irtkImageFunction.h>

#include <irtkGaussianInterpolateImageFunction2D.h>

irtkGaussianInterpolateImageFunction2D::irtkGaussianInterpolateImageFunction2D(double Sigma)
{
  _Sigma = Sigma;
}

irtkGaussianInterpolateImageFunction2D::~irtkGaussianInterpolateImageFunction2D(void)
{}

const char *irtkGaussianInterpolateImageFunction2D::NameOfClass()
{
  return "irtkGaussianInterpolateImageFunction2D";
}

void irtkGaussianInterpolateImageFunction2D::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction::Initialize();

  // Compute size of image
  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();
  this->_z = this->_input->GetZ();

  // Compute voxel size
  this->_input->GetPixelSize(&_xsize, &_ysize, &_zsize);

  // Compute filter extent
  _ExtentX = round(ceil(3*this->_Sigma/_xsize));
  _ExtentY = round(ceil(3*this->_Sigma/_ysize));

  // Compute domain on which the Gaussian interpolation is defined
  this->_x1 = _ExtentX;
  this->_y1 = _ExtentY;
  this->_x2 = this->_input->GetX() - _ExtentX;
  this->_y2 = this->_input->GetY() - _ExtentY;
}

double irtkGaussianInterpolateImageFunction2D::EvaluateInside(double, double, double, double)
{
  cerr << "irtkGaussianInterpolateImageFunction2D::EvaluateInside: Not implemented" << endl;
  return 0;
}

double irtkGaussianInterpolateImageFunction2D::Evaluate(double x, double y, double z, double time)
{
  double val, sum;
  int x1, x2, y1, y2;

  // Initialize
  val = 0;
  sum = 0;
  x1 = x - _ExtentX;
  x2 = x + _ExtentX;
  y1 = y - _ExtentY;
  y2 = y + _ExtentY;

  // Create Gaussian kernel
  irtkScalarGaussian gaussian(this->_Sigma/_xsize, this->_Sigma/_ysize, this->_Sigma/_zsize, x, y, z);

  // Check whether boundary checking is necessary
  if ((x1 > 0) && (x2 < this->_input->GetX()) &&
      (y1 > 0) && (y2 < this->_input->GetY()) &&
      (z  > 0) && (z  < this->_input->GetZ())) {

    // If no, do fast convolution
    for (y = y1; y <= y2; y++) {
      for (x = x1; x <= x2; x++) {
        double temp = gaussian.Evaluate(x, y, z);
        val += temp * this->_input->GetAsDouble(x, y, z, time);
        sum += temp;
      }
    }
  } else {

    // If yes, do slow convolution which handles boundaries
    for (y = y1; y <= y2; y++) {
      for (x = x1; x <= x2; x++) {
        if ((x >= 0) && (x < this->_input->GetX()) &&
            (y >= 0) && (y < this->_input->GetY()) &&
            (z >= 0) && (z < this->_input->GetZ())) {
          double temp = gaussian.Evaluate(x, y, z);
          val += temp * this->_input->GetAsDouble(x, y, z, time);
          sum += temp;
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

