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

inline double cspline(double x)
{
  double xi, xii, xiii;

  xi  = fabs(x);
  xii = xi*xi;
  xiii = xii*xi;
  if (xi < 1) {
    return (xiii - 2*xii + 1);
  } else {
    if (xi < 2) {
      return (-xiii + 5*xii - 8*xi + 4);
    } else {
      return 0;
    }
  }
}

irtkCSplineInterpolateImageFunction2D::irtkCSplineInterpolateImageFunction2D()
{}

irtkCSplineInterpolateImageFunction2D::~irtkCSplineInterpolateImageFunction2D(void)
{}

const char *irtkCSplineInterpolateImageFunction2D::NameOfClass()
{
  return "irtkCSplineInterpolateImageFunction";
}

void irtkCSplineInterpolateImageFunction2D::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction::Initialize();

  // Check if image is 2D
  if (this->_input->GetZ() != 1) {
    cerr << "irtkCSplineInterpolateImageFunction2D::Initialize(): ";
    cerr << "Input image is not 2D" << endl;
    exit(1);
  }

  // Compute size of image
  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();

  // Compute domain on which the cubic spline is defined
  this->_x1 = 1;
  this->_y1 = 1;
  this->_x2 = this->_input->GetX()-2;
  this->_y2 = this->_input->GetY()-2;

  // Compute min and max values
  this->_input->GetMinMaxAsDouble(&this->_min, &this->_max);
}

double irtkCSplineInterpolateImageFunction2D::Evaluate(double x, double y, double z, double time)
{
  double wx, wy, val, sum;
  int i, j, k, l, m, t;

  i = (int)floor(x);
  j = (int)floor(y);
  k = round(z);
  t = round(time);

  val = 0;
  sum = 0;
  for (l = -1; l < 3; l++) {
    if ((l+i >= 0) && (l+i < this->_x)) {
      wx = cspline(l+i-x);
      for (m = -1; m < 3; m++) {
        if ((m+j >= 0) && (m+j < this->_y)) {
          wy = cspline(m+j-y);
          val += wx*wy*this->_input->GetAsDouble(l+i, m+j, k, t);
          sum += wx*wy;
        }
      }
    }
  }

  if (sum != 0) {
    val /= sum;
  } else {
    val = 0;
  }

  if (this->_clamped) {
    if (val > this->_max) return this->_max;
    if (val < this->_min) return this->_min;
  }

  return(val);
}

double irtkCSplineInterpolateImageFunction2D::EvaluateInside(double x, double y, double z, double time)
{
  double wx, wy, val, sum;
  int i, j, k, l, m, t;

  i = (int)floor(x);
  j = (int)floor(y);
  k = round(z);
  t = round(time);

  val = 0;
  sum = 0;
  for (l = -1; l < 3; l++) {

    wx = cspline(l+i-x);
    for (m = -1; m < 3; m++) {
      wy = cspline(m+j-y);
      val += wx*wy*this->_input->GetAsDouble(l+i, m+j, k, t);
      sum += wx*wy;
    }
  }

  if (sum != 0) {
    val /= sum;
  } else {
    val = 0;
  }

  if (this->_clamped) {
    if (val > this->_max) return this->_max;
    if (val < this->_min) return this->_min;
  }

  return(val);
}
