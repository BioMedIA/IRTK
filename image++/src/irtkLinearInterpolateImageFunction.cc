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

#include <irtkLinearInterpolateImageFunction.h>

irtkLinearInterpolateImageFunction::irtkLinearInterpolateImageFunction()
{}

irtkLinearInterpolateImageFunction::~irtkLinearInterpolateImageFunction(void)
{}

const char *irtkLinearInterpolateImageFunction::NameOfClass()
{
  return "irtkLinearInterpolateImageFunction";
}

void irtkLinearInterpolateImageFunction::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction::Initialize();

  // Compute size of image
  this->_x = this->_input->GetX();
  this->_y = this->_input->GetY();
  this->_z = this->_input->GetZ();

  // Compute domain on which the linear interpolation is defined
  this->_x1 = 0;
  this->_y1 = 0;
  this->_z1 = 0;
  this->_x2 = this->_input->GetX() - 1;
  this->_y2 = this->_input->GetY() - 1;
  this->_z2 = this->_input->GetZ() - 1;

  // Calculate offsets for fast pixel access
  this->_offset1 = 0;
  this->_offset2 = 1;
  this->_offset3 = this->_input->GetX();
  this->_offset4 = this->_input->GetX()+1;
  this->_offset5 = this->_input->GetX()*this->_input->GetY();
  this->_offset6 = this->_input->GetX()*this->_input->GetY()+1;
  this->_offset7 = this->_input->GetX()*this->_input->GetY()+this->_input->GetX();
  this->_offset8 = this->_input->GetX()*this->_input->GetY()+this->_input->GetX()+1;
}

double irtkLinearInterpolateImageFunction::EvaluateInside(double x, double y, double z, double time)
{
  int i, j, k;
  double t1, t2, u1, u2, v1, v2;

  // Calculated integer coordinates
  i  = int(x);
  j  = int(y);
  k  = int(z);

  // Calculated fractional coordinates
  t1 = x - i;
  u1 = y - j;
  v1 = z - k;
  t2 = 1 - t1;
  u2 = 1 - u1;
  v2 = 1 - v1;

  // Get pointer to data
  switch (this->_input->GetScalarType()) {
  case IRTK_VOXEL_UNSIGNED_SHORT: {
      unsigned short *ptr = (unsigned short *)this->_input->GetScalarPointer(i, j, k, round(time));

      // Linear interpolation
      return (t1 * (u2 * (v2 * ptr[this->_offset2] + v1 * ptr[this->_offset6]) +
                    u1 * (v2 * ptr[this->_offset4] + v1 * ptr[this->_offset8])) +
              t2 * (u2 * (v2 * ptr[this->_offset1] + v1 * ptr[this->_offset5]) +
                    u1 * (v2 * ptr[this->_offset3] + v1 * ptr[this->_offset7])));

      break;
    }
  case IRTK_VOXEL_SHORT: {
      short *ptr = (short *)this->_input->GetScalarPointer(i, j, k, round(time));

      // Linear interpolation
      return (t1 * (u2 * (v2 * ptr[this->_offset2] + v1 * ptr[this->_offset6]) +
                    u1 * (v2 * ptr[this->_offset4] + v1 * ptr[this->_offset8])) +
              t2 * (u2 * (v2 * ptr[this->_offset1] + v1 * ptr[this->_offset5]) +
                    u1 * (v2 * ptr[this->_offset3] + v1 * ptr[this->_offset7])));

      break;
    }
  case IRTK_VOXEL_FLOAT: {
      float *ptr = (float *)this->_input->GetScalarPointer(i, j, k, round(time));

      // Linear interpolation
      return (t1 * (u2 * (v2 * ptr[this->_offset2] + v1 * ptr[this->_offset6]) +
                    u1 * (v2 * ptr[this->_offset4] + v1 * ptr[this->_offset8])) +
              t2 * (u2 * (v2 * ptr[this->_offset1] + v1 * ptr[this->_offset5]) +
                    u1 * (v2 * ptr[this->_offset3] + v1 * ptr[this->_offset7])));

      break;
    }
  case IRTK_VOXEL_DOUBLE: {
      double *ptr = (double *)this->_input->GetScalarPointer(i, j, k, round(time));

      // Linear interpolation
      return (t1 * (u2 * (v2 * ptr[this->_offset2] + v1 * ptr[this->_offset6]) +
                    u1 * (v2 * ptr[this->_offset4] + v1 * ptr[this->_offset8])) +
              t2 * (u2 * (v2 * ptr[this->_offset1] + v1 * ptr[this->_offset5]) +
                    u1 * (v2 * ptr[this->_offset3] + v1 * ptr[this->_offset7])));

      break;
    }
  default:
    cerr << "irtkLinearInterpolateImageFunction::EvaluateInside: Unknown scalar type" << endl;
    exit(1);
  }
}

double irtkLinearInterpolateImageFunction::Evaluate(double x, double y, double z, double time)
{
  double val;
  double weight;
  int i, j, k, l, m, n, t;

  i = (int)floor(x);
  j = (int)floor(y);
  k = (int)floor(z);
  t = round(time);

  val = 0;
  weight = 0;
  for (l = i; l <= i+1; l++) {
    if ((l >= 0) && (l < this->_x)) {
      for (m = j; m <= j+1; m++) {
        if ((m >= 0) && (m < this->_y)) {
          for (n = k; n <= k+1; n++) {
            if ((n >= 0) && (n < this->_z)) {
              val += (1 - fabs(l - x))*(1 - fabs(m - y))*(1 - fabs(n - z))*this->_input->GetAsDouble(l, m, n, t);
              weight += (1 - fabs(l - x))*(1 - fabs(m - y))*(1 - fabs(n - z));
            }
          }
        }
      }
    }
  }
  
  if (weight > 0)
    val /= weight;
  
  switch (this->_input->GetScalarType()) {
	  case IRTK_VOXEL_UNSIGNED_SHORT: {
		  val = round(val);
		  break;
									  }
	  case IRTK_VOXEL_SHORT: {
		  val = round(val);
		  break;
							 }
	  default: break;
  }
  return val;
}
