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

#include <irtkGaussianBlurring.h>

#include <irtkConvolution.h>

#include <irtkScalarFunctionToImage.h>

template <class VoxelType> irtkGaussianBlurringWithPadding<VoxelType>::irtkGaussianBlurringWithPadding(double Sigma, VoxelType PaddingValue) : irtkGaussianBlurring<VoxelType>(Sigma)
{
  _PaddingValue = PaddingValue;
}

template <class VoxelType> bool irtkGaussianBlurringWithPadding<VoxelType>::RequiresBuffering(void)
{
  return false;
}

template <class VoxelType> const char *irtkGaussianBlurringWithPadding<VoxelType>::NameOfClass()
{
  return "irtkGaussianBlurringWithPadding";
}

template <class VoxelType> void irtkGaussianBlurringWithPadding<VoxelType>::Run()
{
  double xsize, ysize, zsize;

  // Do the initial set up
  this->Initialize();

  // Get voxel dimensions
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);

  // Create scalar function which corresponds to a 1D Gaussian function in X
  irtkScalarGaussian gaussianX(this->_Sigma/xsize, 1, 1, 0, 0, 0);

  // Create filter kernel for 1D Gaussian function in X
  irtkGenericImage<irtkRealPixel> kernelX(2*round(4*this->_Sigma/xsize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceX;
  gaussianSourceX.SetInput (&gaussianX);
  gaussianSourceX.SetOutput(&kernelX);
  gaussianSourceX.Run();

  // Do convolution
  irtkConvolutionWithPadding_1D<VoxelType> convolutionX(this->_PaddingValue);
  convolutionX.SetInput ( this->_input);
  convolutionX.SetInput2(&kernelX);
  convolutionX.SetOutput(this->_output);
  convolutionX.SetNormalization(true);
  convolutionX.irtkImageToImage<VoxelType>::Run();

  // Flip x and y axis of image
  this->_output->FlipXY(1);

  // Create scalar function which corresponds to a 1D Gaussian function in Y
  irtkScalarGaussian gaussianY(this->_Sigma/ysize, 1, 1, 0, 0, 0);

  // Create filter kernel for 1D Gaussian function in Y
  irtkGenericImage<irtkRealPixel> kernelY(2*round(4*this->_Sigma/ysize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceY;
  gaussianSourceY.SetInput (&gaussianY);
  gaussianSourceY.SetOutput(&kernelY);
  gaussianSourceY.Run();

  // Do convolution
  irtkConvolutionWithPadding_1D<VoxelType> convolutionY(this->_PaddingValue);
  convolutionY.SetInput (this->_output);
  convolutionY.SetInput2(&kernelY);
  convolutionY.SetOutput(this->_output);
  convolutionY.SetNormalization(true);
  convolutionY.irtkImageToImage<VoxelType>::Run();

  // Flip x and z axis of image
  this->_output->FlipXZ(1);

  if (this->_output->GetX() != 1) {
    // Create scalar function which corresponds to a 1D Gaussian function in Z
    irtkScalarGaussian gaussianZ(this->_Sigma/zsize, 1, 1, 0, 0, 0);

    // Create filter kernel for 1D Gaussian function in Z
    irtkGenericImage<irtkRealPixel> kernelZ(2*round(4*this->_Sigma/zsize)+1, 1, 1);

    // Do conversion from  scalar function to filter kernel
    irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceZ;
    gaussianSourceZ.SetInput (&gaussianZ);
    gaussianSourceZ.SetOutput(&kernelZ);
    gaussianSourceZ.Run();

    // Do convolution
    irtkConvolutionWithPadding_1D<VoxelType> convolutionZ(this->_PaddingValue);
    convolutionZ.SetInput (this->_output);
    convolutionZ.SetInput2(&kernelZ);
    convolutionZ.SetOutput(this->_output);
    convolutionZ.SetNormalization(true);
    convolutionZ.irtkImageToImage<VoxelType>::Run();
  }

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ(1);
  this->_output->FlipXY(1);

  // Do the final cleaning up
  this->Finalize();
}

template class irtkGaussianBlurringWithPadding<unsigned char>;
template class irtkGaussianBlurringWithPadding<short>;
template class irtkGaussianBlurringWithPadding<unsigned short>;
template class irtkGaussianBlurringWithPadding<float>;
template class irtkGaussianBlurringWithPadding<double>;
