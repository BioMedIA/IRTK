/*=========================================================================

  Date: 14.09.2009
  Author: Laurent Risser
  +++ Works very well +++

=========================================================================*/


#include <irtkImage.h>

///Fast Fourier Transform of numerical recipies (slighly modified)
void four1NR(float data[], unsigned long nn, int isign);

/// * Fast Fourier transform of the complex image contained in 'RealSignal' and 'ImaginarySignal'
///(real and imaginary part).
// * RealSignal and ImaginarySignal MUST have the same size and the size on each dimension MUST
// be a power of 2. (NOT TESTED IN THE FUNCTION)
// * Remark: The FFT of Numerical recipies is called here. Using this function the treated
// signal at the identifier "i" must be in "2*i+1" for the real part and (2*i+2) for the
// imaginary part. If the size of the treated signal is 'S' the size of the input vector
// is "2*s+1" and its value at '0' is not treated (strange notations of numerical recipies!!!).
void DirectFFT(irtkGenericImage<float> * RealSignal,irtkGenericImage<float> * ImaginarySignal);

/// * Inverse Fast Fourier transform of the complex image contained in 'RealSignal' and 'ImaginarySignal'
/// (real and imaginary part).
// * RealSignal and ImaginarySignal MUST have the same size and the size on each dimension MUST
// be a power of 2. (NOT TESTED IN THE FUNCTION)
// * Remark: The FFT of Numerical recipies is called here. Using this function the treated
// signal at the identifier "i" must be in "2*i+1" for the real part and (2*i+2) for the
// imaginary part. If the size of the treated signal is 'S' the size of the input vector
// is "2*s+1" and its value at '0' is not treated (strange notations of numerical recipies!!!).
void InverseFFT(irtkGenericImage<float> * RealSignal,irtkGenericImage<float> * ImaginarySignal);



/// * Convolution in Fourier spaces of the 3D complex image in ('RealPartSignal','ImaginaryPartSignal')
/// by the complex filter in ('RealPartFilter','ImaginaryPartFilter')
// * If the image and/or filter are real, set all imaginary values to 0.
// * The center of the filter is at the coordinate (0,0,0) AND the filter is periodic AND the sum of
//   its values must be 1.
//    => An example of centered filter is then
//     RealPartFilter->Put(0,0,0,0,1./7.);      ImaginaryPartFilter->Put(0,0,0,0,0.);
//     RealPartFilter->Put(1,0,0,0,1./7.);      ImaginaryPartFilter->Put(1,0,0,0,0.);
//     RealPartFilter->Put(0,1,0,0,1./7.);      ImaginaryPartFilter->Put(0,1,0,0,0.);
//     RealPartFilter->Put(0,0,1,0,1./7.);      ImaginaryPartFilter->Put(0,0,1,0,0.);
//     RealPartFilter->Put(NBX-1,0,0,0,1./7.);  ImaginaryPartFilter->Put(NBX-1,0,0,0,0.);
//     RealPartFilter->Put(0,NBY-1,0,0,1./7.);  ImaginaryPartFilter->Put(0,NBY-1,0,0,0.);
//     RealPartFilter->Put(0,0,NBZ-1,0,1./7.);  ImaginaryPartFilter->Put(0,0,NBZ-1,0,0.);
// * RealPartSignal, ImaginaryPartSignal, RealPartFilter and ImaginaryPartFilter MUST have the same size 
// and the size on each dimension MUST be a power of 2. (NOT TESTED IN THE FUNCTION)
void ConvolutionInFourier(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter);


/// * Convolution in Fourier spaces of the 3D complex image in ('RealPartSignal','ImaginaryPartSignal')
/// by complex filter in ('RealPartFilter','ImaginaryPartFilter') which is ALREADY in Fourier space
void ConvolutionInFourierNoFilterTransfo(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilterTransformedFrSpace,irtkGenericImage<float> * ImaginaryPartFilterTransformedFrSpace);






/// * Deconvolution in Fourier spaces of the 3D complex image in ('RealPartSignal','ImaginaryPartSignal')
/// by the complex filter in ('RealPartFilter','ImaginaryPartFilter')
// * RealPartSignal, ImaginaryPartSignal, RealPartFilter and ImaginaryPartFilter MUST have the same size 
// and the size on each dimension MUST be a power of 2. (NOT TESTED IN THE FUNCTION)
void DeconvolutionInFourier(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter);


/// * Deconvolution in Fourier spaces of the 3D complex image in ('RealPartSignal','ImaginaryPartSignal')
/// by complex filter in ('RealPartFilter','ImaginaryPartFilter') which is ALREADY in Fourier space
void DeconvolutionInFourierNoFilterTransfo(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilterTransformedFrSpace,irtkGenericImage<float> * ImaginaryPartFilterTransformedFrSpace);


/// Build a Gaussian filter of standard deviation 'sigma'
void MakeGaussianFilter(float sigma,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter);

/// Build a Gaussian filter of aniotropic standard deviation 'sigmaX,sigmaY,sigmaZ'
void MakeAnisotropicGaussianFilter(float weight,float sigmaX,float sigmaY,float sigmaZ,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter);

/// Sum of aniotropic Gaussian filters
void MakeSumOf2AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter);

void MakeSumOf3AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter);

void MakeSumOf4AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter);
