%include <irtkVoxel.h>

%include "irtkImageAttributes.i"
%include "irtkBaseImage.i"
%include "irtkGenericImage.i"
%include "irtkGaussianBlurring.i"

%template(irtkByteImage) irtkGenericImage<irtkBytePixel>;
%template(irtkGreyImage) irtkGenericImage<irtkGreyPixel>;
%template(irtkRealImage) irtkGenericImage<irtkRealPixel>;
//// %template(irtkVector3DCharImage)  irtkGenericImage<irtkVector3D<char> >;
//// %template(irtkVector3DShortImage)  irtkGenericImage<irtkVector3D<short> >;
//// %template(irtkVector3DFloatImage)  irtkGenericImage<irtkVector3D<float> >;
//// %template(irtkVector3DDoubleImage)  irtkGenericImage<irtkVector3D<double> >;
typedef irtkGenericImage<irtkBytePixel> irtkByteImage;
typedef irtkGenericImage<irtkGreyPixel> irtkGreyImage;
typedef irtkGenericImage<irtkRealPixel> irtkRealImage;
//// typedef irtkGenericImage<irtkVector3D<char> >  irtkVector3DImage;
//// typedef irtkGenericImage<irtkVector3D<short> > irtkVector3DShortImage;
//// typedef irtkGenericImage<irtkVector3D<float> > irtkVector3DFloatImage;
//// typedef irtkGenericImage<irtkVector3D<double> > irtkVector3DDoubleImage;

%template(irtkByteGaussianBlurring) irtkGaussianBlurring<irtkBytePixel>;
%template(irtkGreyGaussianBlurring) irtkGaussianBlurring<irtkGreyPixel>;
%template(irtkRealGaussianBlurring) irtkGaussianBlurring<irtkRealPixel>;
typedef irtkGaussianBlurring<irtkBytePixel> irtkByteGaussianBlurring;
typedef irtkGaussianBlurring<irtkGreyPixel> irtkGreyGaussianBlurring;
typedef irtkGaussianBlurring<irtkRealPixel> irtkRealGaussianBlurring;
