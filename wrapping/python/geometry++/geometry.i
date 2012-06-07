%include "irtkPoint.i"
%include "irtkMatrix.i"
%include "irtkVector.i"
%include "irtkVector3D.i"

%template(irtkVector3DChar) irtkVector3D<char>;
%template(irtkVector3DShort) irtkVector3D<short>;
%template(irtkVector3DFloat) irtkVector3D<float>;
%template(irtkVector3DDouble) irtkVector3D<double>;
typedef irtkVector3D<char> irtkVector3DChar;
typedef irtkVector3D<short> irtkVector3DShort;
typedef irtkVector3D<float> irtkVector3DFloat;
typedef irtkVector3D<double> irtkVector3DDouble;
