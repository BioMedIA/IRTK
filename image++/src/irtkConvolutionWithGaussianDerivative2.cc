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

#include <irtkConvolution.h>
#include <irtkScalarFunctionToImage.h>
#include <irtkScalarGaussianDx.h>
#include <irtkScalarGaussianDy.h>
#include <irtkScalarGaussianDz.h>
#include <irtkScalarGaussianDxDx.h>
#include <irtkScalarGaussianDyDy.h>
#include <irtkScalarGaussianDzDz.h>
#include <irtkScalarGaussianDxDz.h>
#include <irtkScalarGaussianDyDz.h>
#include <irtkScalarGaussianDxDy.h>

#include <irtkConvolutionWithGaussianDerivative2.h>

template <class VoxelType> irtkConvolutionWithGaussianDerivative2<VoxelType>::irtkConvolutionWithGaussianDerivative2(double Sigma)
{
  _Sigma = Sigma;
}

template <class VoxelType> irtkConvolutionWithGaussianDerivative2<VoxelType>::~irtkConvolutionWithGaussianDerivative2()
{}

template <class VoxelType> const char *irtkConvolutionWithGaussianDerivative2<VoxelType>::NameOfClass()
{
  return "irtkConvolutionWithGaussianDerivative2";
}

template <class VoxelType> bool irtkConvolutionWithGaussianDerivative2<VoxelType>::RequiresBuffering()
{
  return true;
}

template <class VoxelType> void irtkConvolutionWithGaussianDerivative2<VoxelType>::Ixx()
{

  double xsize, ysize, zsize;


  // Do the initial set up
  this->Initialize();

  // Get voxel dimensions
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);

  // Create scalar function which corresponds to a 1D Gaussian function in X
  irtkScalarGaussianDxDx gaussianDxDx(this->_Sigma/xsize, 1, 1, 0, 0, 0);
  // Create filter kernel for 1D Gaussian function in X
  irtkGenericImage<irtkRealPixel> kernelX(2*round(4*this->_Sigma/xsize)+1, 1, 1);
  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceX;
  gaussianSourceX.SetInput (&gaussianDxDx);
  gaussianSourceX.SetOutput(&kernelX);
  gaussianSourceX.Run();

  // Do convolution
  irtkConvolution_1D<VoxelType> convolutionX;
  convolutionX.SetInput ( this->_input);
  convolutionX.SetInput2(&kernelX);
  convolutionX.SetOutput(this->_output);
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

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ();
  this->_output->FlipXY();

  // Do the final cleaning up
  this->Finalize();

}

template <class VoxelType> void irtkConvolutionWithGaussianDerivative2<VoxelType>::Ixy()
{
  double xsize, ysize, zsize;

  // Do the initial set up
  this->Initialize();

  // Get voxel dimensions
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);

  // Create scalar function which corresponds to a 1D Gaussian function in X
  irtkScalarGaussianDx gaussianDx(this->_Sigma/xsize, 1, 1, 0, 0, 0);
  // Create filter kernel for 1D Gaussian function in X
  irtkGenericImage<irtkRealPixel> kernelX(2*round(4*this->_Sigma/xsize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceX;
  gaussianSourceX.SetInput (&gaussianDx);
  gaussianSourceX.SetOutput(&kernelX);
  gaussianSourceX.Run();

  // Do convolution
  irtkConvolution_1D<VoxelType> convolutionX;
  convolutionX.SetInput ( this->_input);
  convolutionX.SetInput2(&kernelX);
  convolutionX.SetOutput(this->_output);
  convolutionX.irtkImageToImage<VoxelType>::Run();

  // Flip x and y axis of image
  this->_output->FlipXY();

  // Create scalar function which corresponds to a 1D Gaussian function in Y
  irtkScalarGaussianDx gaussianDy(this->_Sigma/ysize, 1, 1, 0, 0, 0);

  // Create filter kernel for 1D Gaussian function in Y
  irtkGenericImage<irtkRealPixel> kernelY(2*round(4*this->_Sigma/ysize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceY;
  gaussianSourceY.SetInput (&gaussianDy);
  gaussianSourceY.SetOutput(&kernelY);
  gaussianSourceY.Run();

  // Do convolution
  irtkConvolution_1D<VoxelType> convolutionY;
  convolutionY.SetInput (this->_output);
  convolutionY.SetInput2(&kernelY);
  convolutionY.SetOutput(this->_output);
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

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ();
  this->_output->FlipXY();

  // Do the final cleaning up
  this->Finalize();

}

template <class VoxelType> void irtkConvolutionWithGaussianDerivative2<VoxelType>::Ixz()
{
  double xsize, ysize, zsize;

  // Do the initial set up
  this->Initialize();

  // Get voxel dimensions
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);

  // Create scalar function which corresponds to a 1D Gaussian function in X
  irtkScalarGaussianDx gaussianDx(this->_Sigma/xsize, 1, 1, 0, 0, 0);
  // Create filter kernel for 1D Gaussian function in X
  irtkGenericImage<irtkRealPixel> kernelX(2*round(4*this->_Sigma/xsize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceX;
  gaussianSourceX.SetInput (&gaussianDx);
  gaussianSourceX.SetOutput(&kernelX);
  gaussianSourceX.Run();

  // Do convolution
  irtkConvolution_1D<VoxelType> convolutionX;
  convolutionX.SetInput ( this->_input);
  convolutionX.SetInput2(&kernelX);
  convolutionX.SetOutput(this->_output);
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
    irtkScalarGaussianDx gaussianDz(this->_Sigma/zsize, 1, 1, 0, 0, 0);

    // Create filter kernel for 1D Gaussian function in Z
    irtkGenericImage<irtkRealPixel> kernelZ(2*round(4*this->_Sigma/zsize)+1, 1, 1);

    // Do conversion from  scalar function to filter kernel
    irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceZ;
    gaussianSourceZ.SetInput (&gaussianDz);
    gaussianSourceZ.SetOutput(&kernelZ);
    gaussianSourceZ.Run();

    // Do convolution
    irtkConvolution_1D<VoxelType> convolutionZ;
    convolutionZ.SetInput (this->_output);
    convolutionZ.SetInput2(&kernelZ);
    convolutionZ.SetOutput(this->_output);
    convolutionZ.irtkImageToImage<VoxelType>::Run();
  }

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ();
  this->_output->FlipXY();

  // Do the final cleaning up
  this->Finalize();
}

template <class VoxelType> void irtkConvolutionWithGaussianDerivative2<VoxelType>::Iyy()
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
  convolutionX.SetNormalization(true);
  convolutionX.irtkImageToImage<VoxelType>::Run();

  // Flip x and y axis of image
  this->_output->FlipXY();

  // Create scalar function which corresponds to a 1D Gaussian function in Y
  irtkScalarGaussianDxDx gaussianDyDy(this->_Sigma/ysize, 1, 1, 0, 0, 0);

  // Create filter kernel for 1D Gaussian function in Y
  irtkGenericImage<irtkRealPixel> kernelY(2*round(4*this->_Sigma/ysize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceY;
  gaussianSourceY.SetInput (&gaussianDyDy);
  gaussianSourceY.SetOutput(&kernelY);
  gaussianSourceY.Run();

  // Do convolution
  irtkConvolution_1D<VoxelType> convolutionY;
  convolutionY.SetInput (this->_output);
  convolutionY.SetInput2(&kernelY);
  convolutionY.SetOutput(this->_output);
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

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ();
  this->_output->FlipXY();

  // Do the final cleaning up
  this->Finalize();
}

template <class VoxelType> void irtkConvolutionWithGaussianDerivative2<VoxelType>::Iyz()
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
  convolutionX.SetNormalization(true);
  convolutionX.irtkImageToImage<VoxelType>::Run();

  // Flip x and y axis of image
  this->_output->FlipXY();

  // Create scalar function which corresponds to a 1D Gaussian function in Y
  irtkScalarGaussianDx gaussianDy(this->_Sigma/ysize, 1, 1, 0, 0, 0);

  // Create filter kernel for 1D Gaussian function in Y
  irtkGenericImage<irtkRealPixel> kernelY(2*round(4*this->_Sigma/ysize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceY;
  gaussianSourceY.SetInput (&gaussianDy);
  gaussianSourceY.SetOutput(&kernelY);
  gaussianSourceY.Run();

  // Do convolution
  irtkConvolution_1D<VoxelType> convolutionY;
  convolutionY.SetInput (this->_output);
  convolutionY.SetInput2(&kernelY);
  convolutionY.SetOutput(this->_output);
  convolutionY.irtkImageToImage<VoxelType>::Run();

  // Flip x and z axis of image
  this->_output->FlipXZ();

  if (this->_output->GetZ() != 1) {
    // Create scalar function which corresponds to a 1D Gaussian function in Z
    irtkScalarGaussianDx gaussianDz(this->_Sigma/zsize, 1, 1, 0, 0, 0);

    // Create filter kernel for 1D Gaussian function in Z
    irtkGenericImage<irtkRealPixel> kernelZ(2*round(4*this->_Sigma/zsize)+1, 1, 1);

    // Do conversion from  scalar function to filter kernel
    irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceZ;
    gaussianSourceZ.SetInput (&gaussianDz);
    gaussianSourceZ.SetOutput(&kernelZ);
    gaussianSourceZ.Run();

    // Do convolution
    irtkConvolution_1D<VoxelType> convolutionZ;
    convolutionZ.SetInput (this->_output);
    convolutionZ.SetInput2(&kernelZ);
    convolutionZ.SetOutput(this->_output);
    convolutionZ.irtkImageToImage<VoxelType>::Run();
  }

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ();
  this->_output->FlipXY();

  // Do the final cleaning up
  this->Finalize();
}

template <class VoxelType> void irtkConvolutionWithGaussianDerivative2<VoxelType>::Izz()
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
    irtkScalarGaussianDxDx gaussianDzDz(this->_Sigma/zsize, 1, 1, 0, 0, 0);

    // Create filter kernel for 1D Gaussian function in Z
    irtkGenericImage<irtkRealPixel> kernelZ(2*round(4*this->_Sigma/zsize)+1, 1, 1);

    // Do conversion from  scalar function to filter kernel
    irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceZ;
    gaussianSourceZ.SetInput (&gaussianDzDz);
    gaussianSourceZ.SetOutput(&kernelZ);
    gaussianSourceZ.Run();

    // Do convolution
    irtkConvolution_1D<VoxelType> convolutionZ;
    convolutionZ.SetInput (this->_output);
    convolutionZ.SetInput2(&kernelZ);
    convolutionZ.SetOutput(this->_output);
    convolutionZ.irtkImageToImage<VoxelType>::Run();
  }

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ();
  this->_output->FlipXY();

  // Do the final cleaning up
  this->Finalize();
}


template class irtkConvolutionWithGaussianDerivative2<irtkBytePixel>;
template class irtkConvolutionWithGaussianDerivative2<irtkGreyPixel>;
template class irtkConvolutionWithGaussianDerivative2<irtkRealPixel>;
