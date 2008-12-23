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

template <class VoxelType> irtkGaussianBlurring<VoxelType>::irtkGaussianBlurring(double Sigma)
{
  _Sigma = Sigma;
}

template <class VoxelType> irtkGaussianBlurring<VoxelType>::~irtkGaussianBlurring(void)
{
}

template <class VoxelType> Bool irtkGaussianBlurring<VoxelType>::RequiresBuffering(void)
{
  return False;
}

template <class VoxelType> const char *irtkGaussianBlurring<VoxelType>::NameOfClass()
{
  return "irtkGaussianBlurring";
}

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::Run()
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
  irtkConvolution_1D<VoxelType> convolutionX;
  convolutionX.SetInput ( this->_input);
  convolutionX.SetInput2(&kernelX);
  convolutionX.SetOutput(this->_output);
  convolutionX.SetNormalization(True);
  convolutionX.irtkImageToImage<VoxelType>::Run();

  // Flip x and y axis of image
  this->_output->FlipXY();

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
  irtkConvolution_1D<VoxelType> convolutionY;
  convolutionY.SetInput (this->_output);
  convolutionY.SetInput2(&kernelY);
  convolutionY.SetOutput(this->_output);
  convolutionY.SetNormalization(True);
  convolutionY.irtkImageToImage<VoxelType>::Run();

  // Flip x and z axis of image
  this->_output->FlipXZ();

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
    irtkConvolution_1D<VoxelType> convolutionZ;
    convolutionZ.SetInput (this->_output);
    convolutionZ.SetInput2(&kernelZ);
    convolutionZ.SetOutput(this->_output);
    convolutionZ.SetNormalization(True);
    convolutionZ.irtkImageToImage<VoxelType>::Run();
  }

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ();
  this->_output->FlipXY();

  // Do the final cleaning up
  this->Finalize();
}

template class irtkGaussianBlurring<irtkBytePixel>;
template class irtkGaussianBlurring<irtkGreyPixel>;
template class irtkGaussianBlurring<irtkRealPixel>;
