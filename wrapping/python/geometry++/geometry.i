%include "itkPoint.i"
%include "itkMatrix.i"
%include "itkQuaternion.i"
%include "itkVector.i"
%include "itkVector3D.i"

%template(itkVector3DChar) itkVector3D<char>;
%template(itkVector3DShort) itkVector3D<short>;
%template(itkVector3DFloat) itkVector3D<float>;
%template(itkVector3DDouble) itkVector3D<double>;
typedef itkVector3D<char> itkVector3DChar;
typedef itkVector3D<short> itkVector3DShort;
typedef itkVector3D<float> itkVector3DFloat;
typedef itkVector3D<double> itkVector3DDouble;
