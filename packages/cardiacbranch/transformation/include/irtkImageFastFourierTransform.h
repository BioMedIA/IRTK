/*=========================================================================

  Date: 14.09.2009
  Author: Laurent Risser
  +++ Works very well +++

=========================================================================*/


#include <irtkImage.h>

//Fast Fourier Transform of numerical recipies (slighly modified)
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




/// * Deconvolution in Fourier spaces of the 3D complex image in ('RealPartSignal','ImaginaryPartSignal')
/// by the complex filter in ('RealPartFilter','ImaginaryPartFilter')
// * RealPartSignal, ImaginaryPartSignal, RealPartFilter and ImaginaryPartFilter MUST have the same size 
// and the size on each dimension MUST be a power of 2. (NOT TESTED IN THE FUNCTION)
void DeconvolutionInFourier(irtkGenericImage<float> * RealPartSignal,irtkGenericImage<float> * ImaginaryPartSignal,irtkGenericImage<float> * RealPartFilter,irtkGenericImage<float> * ImaginaryPartFilter);
