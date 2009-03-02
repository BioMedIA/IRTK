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

irtkNearestNeighborInterpolateImageFunction2D::irtkNearestNeighborInterpolateImageFunction2D()
{
}

irtkNearestNeighborInterpolateImageFunction2D::~irtkNearestNeighborInterpolateImageFunction2D(void)
{
}

const char *irtkNearestNeighborInterpolateImageFunction2D::NameOfClass()
{
  return "irtkNearestNeighborInterpolateImageFunction2D";
}

void irtkNearestNeighborInterpolateImageFunction2D::Initialize()
{
  /// Initialize baseclass
  this->irtkImageFunction::Initialize();

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
  this->_x1 = -0.5;
  this->_y1 = -0.5;
  this->_x2 = this->_input->GetX()-0.5;
  this->_y2 = this->_input->GetY()-0.5;
}

double irtkNearestNeighborInterpolateImageFunction2D::EvaluateInside(double x, double y, double z, double t)
{
  int i, j, k, l;

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  return this->_input->GetAsDouble(i, j, k, l);
}

double irtkNearestNeighborInterpolateImageFunction2D::Evaluate(double x, double y, double z, double t)
{
  int i, j, k, l;

  i = round(x);
  j = round(y);
  k = round(z);
  l = round(t);

  if ((i < 0) || (i >= this->_x) || (j < 0) || (j >= this->_y)) {
    return this->_DefaultValue;
  } else {
    return this->_input->GetAsDouble(i, j, k, l);
  }
}

