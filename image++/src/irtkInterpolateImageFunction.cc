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

irtkInterpolateImageFunction::irtkInterpolateImageFunction() : irtkImageFunction()
{
  // Default parameters
  _clamped = true;
}

irtkInterpolateImageFunction::~irtkInterpolateImageFunction()
{
  // Set input
  this->_input  = NULL;
}

irtkInterpolateImageFunction *irtkInterpolateImageFunction::New(irtkInterpolationMode interpolationMode, irtkBaseImage *image)
{
  if (image->GetZ() == 1) {
    switch (interpolationMode) {
    case Interpolation_NN:
      return new irtkNearestNeighborInterpolateImageFunction2D;
      break;
    case Interpolation_Linear:
      return new irtkLinearInterpolateImageFunction2D;
      break;
    case Interpolation_CSpline:
      return new irtkCSplineInterpolateImageFunction2D;
      break;
    case Interpolation_BSpline:
      return new irtkBSplineInterpolateImageFunction2D;
      break;
    case Interpolation_Sinc:
      return new irtkSincInterpolateImageFunction2D;
      break;
    case Interpolation_Gaussian:
      return new irtkGaussianInterpolateImageFunction2D;
      break;
    }
  } else {
    switch (interpolationMode) {
    case Interpolation_NN:
      return new irtkNearestNeighborInterpolateImageFunction;
      break;
    case Interpolation_Linear:
      return new irtkLinearInterpolateImageFunction;
      break;
    case Interpolation_CSpline:
      return new irtkCSplineInterpolateImageFunction;
      break;
    case Interpolation_BSpline:
      return new irtkBSplineInterpolateImageFunction;
      break;
    case Interpolation_Sinc:
      return new irtkSincInterpolateImageFunction;
      break;
    case Interpolation_Gaussian:
      return new irtkGaussianInterpolateImageFunction;
      break;
    default:
      break;
    }
  }
  return NULL;
}

void irtkInterpolateImageFunction::Initialize()
{
  // Initialize base class
  this->irtkImageFunction::Initialize();

  // Set default input image domain
  this->_z1 = -0.5;
  this->_z2 = +0.5;
}

