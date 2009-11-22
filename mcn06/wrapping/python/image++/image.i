%include <itkBytePixel.h>
%include <itkGreyPixel.h>
%include <itkRealPixel.h>

%include "itkBaseImage.i"
%include "itkGenericImage.i"
%include "itkGaussianBlurring.i"

%template(itkByteImage) itkGenericImage<itkBytePixel>;
%template(itkGreyImage) itkGenericImage<itkGreyPixel>;
%template(itkRealImage) itkGenericImage<itkRealPixel>;
%template(itkVector3DCharImage)  itkGenericImage<itkVector3D<char> >;
%template(itkVector3DShortImage)  itkGenericImage<itkVector3D<short> >;
%template(itkVector3DFloatImage)  itkGenericImage<itkVector3D<float> >;
%template(itkVector3DDoubleImage)  itkGenericImage<itkVector3D<double> >;
typedef itkGenericImage<itkBytePixel> itkByteImage;
typedef itkGenericImage<itkGreyPixel> itkGreyImage;
typedef itkGenericImage<itkRealPixel> itkRealImage;
typedef itkGenericImage<itkVector3D<char> >  itkVector3DImage;
typedef itkGenericImage<itkVector3D<short> > itkVector3DShortImage;
typedef itkGenericImage<itkVector3D<float> > itkVector3DFloatImage;
typedef itkGenericImage<itkVector3D<double> > itkVector3DDoubleImage;

%template(itkByteGaussianBlurring) itkGaussianBlurring<itkBytePixel>;
%template(itkGreyGaussianBlurring) itkGaussianBlurring<itkGreyPixel>;
%template(itkRealGaussianBlurring) itkGaussianBlurring<itkRealPixel>;
typedef itkGaussianBlurring<itkBytePixel> itkByteGaussianBlurring;
typedef itkGaussianBlurring<itkGreyPixel> itkGreyGaussianBlurring;
typedef itkGaussianBlurring<itkRealPixel> itkRealGaussianBlurring;
