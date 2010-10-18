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

#include <irtkGaussianBlurring4D.h>

#include <irtkConvolution.h>

#include <irtkScalarFunctionToImage.h>

template <class VoxelType> irtkGaussianBlurring4D<VoxelType>::irtkGaussianBlurring4D(double Sigma)
{
  _Sigma = Sigma;
}

template <class VoxelType> irtkGaussianBlurring4D<VoxelType>::~irtkGaussianBlurring4D(void)
{
}

template <class VoxelType> bool irtkGaussianBlurring4D<VoxelType>::RequiresBuffering(void)
{
  return false;
}

template <class VoxelType> const char *irtkGaussianBlurring4D<VoxelType>::NameOfClass()
{
  return "irtkGaussianBlurring4D";
}

template <class VoxelType> void irtkGaussianBlurring4D<VoxelType>::Run()
{
  double xsize, ysize, zsize, tsize;

  // Do the initial set up
  this->Initialize();

  // Get voxel dimensions
  this->_input->GetPixelSize(&xsize, &ysize, &zsize, &tsize);

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
  convolutionX.SetNormalization(true);
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
  convolutionY.SetNormalization(true);
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
    convolutionZ.SetNormalization(true);
    convolutionZ.irtkImageToImage<VoxelType>::Run();
  }

  // Flip x and t axis of image
  this->_output->FlipXT();

  // Create scalar function which corresponds to a 1D Gaussian function in T
  irtkScalarGaussian gaussianT(this->_Sigma/tsize, 1, 1, 0, 0, 0);

  // Create filter kernel for 1D Gaussian function in T
  irtkGenericImage<irtkRealPixel> kernelT(2*round(4*this->_Sigma/tsize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceT;
  gaussianSourceT.SetInput (&gaussianT);
  gaussianSourceT.SetOutput(&kernelT);
  gaussianSourceT.Run();

  // Do convolution
  irtkConvolution_1D<VoxelType> convolutionT;
  convolutionT.SetInput (this->_output);
  convolutionT.SetInput2(&kernelT);
  convolutionT.SetOutput(this->_output);
  convolutionT.SetNormalization(true);
  convolutionT.irtkImageToImage<VoxelType>::Run();

  // Flip image back, first t and z axis, then x and z axis and finally x and y axis
  this->_output->FlipXT();
  this->_output->FlipXZ();
  this->_output->FlipXY();

  // Do the final cleaning up
  this->Finalize();
}

template class irtkGaussianBlurring4D<unsigned char>;
template class irtkGaussianBlurring4D<short>;
template class irtkGaussianBlurring4D<unsigned short>;
template class irtkGaussianBlurring4D<float>;
template class irtkGaussianBlurring4D<double>;
