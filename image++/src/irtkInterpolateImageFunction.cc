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

template <class VoxelType> irtkInterpolateImageFunction<VoxelType>::irtkInterpolateImageFunction() : irtkImageFunction<VoxelType>()
{
  // Default parameters
  _clamped = True;
}

template <class VoxelType> irtkInterpolateImageFunction<VoxelType>::~irtkInterpolateImageFunction()
{
  // Set input
  this->_input  = NULL;
}

template <class VoxelType> irtkInterpolateImageFunction<VoxelType> *irtkInterpolateImageFunction<VoxelType>::New(irtkInterpolationMode interpolationMode, irtkBaseImage *image)
{
  if (image->GetZ() == 1) {
    switch (interpolationMode) {
    case Interpolation_NN:
      return new irtkNearestNeighborInterpolateImageFunction2D<VoxelType>;
      break;
    case Interpolation_Linear:
      return new irtkLinearInterpolateImageFunction2D<VoxelType>;
      break;
    case Interpolation_CSpline:
      return new irtkCSplineInterpolateImageFunction2D<VoxelType>;
      break;
    case Interpolation_BSpline:
      return new irtkBSplineInterpolateImageFunction2D<VoxelType>;
      break;
    case Interpolation_Sinc:
      return new irtkSincInterpolateImageFunction2D<VoxelType>;
      break;
    case Interpolation_Gaussian:
      return new irtkGaussianInterpolateImageFunction2D<VoxelType>;
      break;
    }
  } else {
    switch (interpolationMode) {
    case Interpolation_NN:
      return new irtkNearestNeighborInterpolateImageFunction<VoxelType>;
      break;
    case Interpolation_Linear:
      return new irtkLinearInterpolateImageFunction<VoxelType>;
      break;
    case Interpolation_CSpline:
      return new irtkCSplineInterpolateImageFunction<VoxelType>;
      break;
    case Interpolation_BSpline:
      return new irtkBSplineInterpolateImageFunction<VoxelType>;
      break;
    case Interpolation_Sinc:
      return new irtkSincInterpolateImageFunction<VoxelType>;
      break;
    case Interpolation_Gaussian:
      return new irtkGaussianInterpolateImageFunction<VoxelType>;
      break;
    default:
      break;
    }
  }
  return NULL;
}

template <class VoxelType> void irtkInterpolateImageFunction<VoxelType>::Initialize()
{
  // Initialize base class
  this->irtkImageFunction<VoxelType>::Initialize();

  // Set default input image domain
  this->_z1 = -0.5;
  this->_z2 = +0.5;
}

template class irtkInterpolateImageFunction<irtkBytePixel>;
template class irtkInterpolateImageFunction<irtkGreyPixel>;
template class irtkInterpolateImageFunction<irtkRealPixel>;
