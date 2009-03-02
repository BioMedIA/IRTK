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

irtkCSplineInterpolateImageFunction::irtkCSplineInterpolateImageFunction()
{}

irtkCSplineInterpolateImageFunction::~irtkCSplineInterpolateImageFunction(void)
{}

const char *irtkCSplineInterpolateImageFunction::NameOfClass()
{
  return "irtkCSplineInterpolateImageFunction";
}

void irtkCSplineInterpolateImageFunction::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction::Initialize();

  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();
  this->_z = this->_input->GetZ();

  // Compute domain on which the cubic spline is defined
  this->_x1 = 1;
  this->_y1 = 1;
  this->_z1 = 1;
  this->_x2 = this->_input->GetX()-2;
  this->_y2 = this->_input->GetY()-2;
  this->_z2 = this->_input->GetZ()-2;

  // Compute min and max values
  this->_input->GetMinMaxAsDouble(&this->_min, &this->_max);
}

double irtkCSplineInterpolateImageFunction::Evaluate(double x, double y, double z, double time)
{
  double wx, wy, wz, val, sum;
  int i, j, k, l, m, n, t;

  i = (int)floor(x);
  j = (int)floor(y);
  k = (int)floor(z);
  t = round(time);

  val = 0;
  sum = 0;
  for (l = -1; l < 3; l++) {
    if ((l+i >= 0) && (l+i < this->_x)) {
      wx = cspline(l+i-x);
      for (m = -1; m < 3; m++) {
        if ((m+j >= 0) && (m+j < this->_y)) {
          wy = cspline(m+j-y);
          for (n = -1; n < 3; n++) {
            if ((n+k >= 0) && (n+k < this->_z)) {
              wz = wx*wy*cspline(n+k-z);
              val += wz*this->_input->GetAsDouble(l+i, m+j, n+k, t);
              sum += wz;
            }
          }
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

double irtkCSplineInterpolateImageFunction::EvaluateInside(double x, double y, double z, double time)
{
  double wx, wy, wz, val, sum;
  int i, j, k, l, m, n, t;

  i = (int)floor(x);
  j = (int)floor(y);
  k = (int)floor(z);
  t = round(time);

  val = 0;
  sum = 0;
  for (l = -1; l < 3; l++) {
    wx = cspline(l+i-x);
    for (m = -1; m < 3; m++) {
      wy = cspline(m+j-y);
      for (n = -1; n < 3; n++) {
        wz = wx*wy*cspline(n+k-z);
        val += wz*this->_input->GetAsDouble(l+i, m+j, n+k, t);
        sum += wz;
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
