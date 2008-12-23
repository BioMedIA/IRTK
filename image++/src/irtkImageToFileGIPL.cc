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

#include <irtkImageToFile.h>

#include <irtkGIPL.h>

template <class VoxelType> const char *irtkImageToFileGIPL<VoxelType>::NameOfClass()
{
  return "irtkImageToFileGIPL";
}

template <> void irtkImageToFileGIPL<irtkBytePixel>::Initialize()
{
  int i;
  irtkBytePixel min, max;
  double xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3], pos[3];

  // Initialize base class
  this->irtkImageToFile<irtkBytePixel>::Initialize();

  // Write image dimensions
  this->WriteAsShort(_input->GetX(), 0);
  this->WriteAsShort(_input->GetY(), 2);
  this->WriteAsShort(_input->GetZ(), 4);
  this->WriteAsShort(1, 6);

  // Write type
  this->WriteAsShort(GIPL_U_CHAR, 8);

  // Write voxel dimensions
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  this->WriteAsFloat(static_cast<float>(xsize), 10);
  this->WriteAsFloat(static_cast<float>(ysize), 14);
  this->WriteAsFloat(static_cast<float>(zsize), 18);
  this->WriteAsFloat(float(0), 22);

  // Write patient description
  for (i = 0; i < 80; i++) {
    this->WriteAsChar(char(0), 26+i*sizeof(char));
  }

  // Write image orientation, checkign for negative determinat first to issue a warnign
  _input->GetOrientation(xaxis, yaxis, zaxis);
  irtkMatrix m(3, 3);
  for (i = 0; i < 3; i++) {
    m(i, 0) = xaxis[i];
    m(i, 1) = yaxis[i];
    m(i, 2) = zaxis[i];
  }
  if (m.Det() < 0) {
    cerr << this->NameOfClass() << "::Initialize() : Unsupported axis orientation in output file\n";
  }
  this->WriteAsFloat(static_cast<float>(xaxis[0]), 106);
  this->WriteAsFloat(static_cast<float>(xaxis[1]), 110);
  this->WriteAsFloat(static_cast<float>(xaxis[2]), 114);
  this->WriteAsFloat(static_cast<float>(yaxis[0]), 122);
  this->WriteAsFloat(static_cast<float>(yaxis[1]), 126);
  this->WriteAsFloat(static_cast<float>(yaxis[2]), 130);

  // Write identifier
  this->WriteAsInt(0, 154);

  // Write spare bytes
  for (i = 0; i < 28; i++) {
    this->WriteAsChar(char(0), 158+i*sizeof(char));
  }

  // Write flag and orientation
  this->WriteAsChar(char(0), 186);
  this->WriteAsChar(char(0), 187);

  // Write min and max values
  _input->GetMinMax(&min, &max);
  this->WriteAsDouble(min, 188);
  this->WriteAsDouble(max, 196);

  // Write origin
  _input->GetOrigin(pos[0], pos[1], pos[2]);
  this->WriteAsDouble(pos[0], 204);
  this->WriteAsDouble(pos[1], 212);
  this->WriteAsDouble(pos[2], 220);
  this->WriteAsDouble(double(0), 228);

  // Write magic number
  this->WriteAsFloat(float(0), 236);
  this->WriteAsFloat(float(0), 240);
  this->WriteAsUInt(GIPL_MAGIC_EXT, 244);
  this->WriteAsFloat(float(0), 248);
  this->WriteAsUInt(GIPL_MAGIC1, 252);

  // Calculate data address
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = GIPL_HEADERSIZE + i * _input->GetX() * _input->GetY() * sizeof(irtkBytePixel);
  }
}

template <> void irtkImageToFileGIPL<irtkGreyPixel>::Initialize()
{
  int i;
  irtkGreyPixel min, max;
  double xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3], pos[3];

  // Initialize base class
  this->irtkImageToFile<irtkGreyPixel>::Initialize();

  // Write image dimensions
  this->WriteAsShort(_input->GetX(), 0);
  this->WriteAsShort(_input->GetY(), 2);
  this->WriteAsShort(_input->GetZ(), 4);
  this->WriteAsShort(1, 6);

  // Write type
  this->WriteAsShort(GIPL_SHORT, 8);

  // Write voxel dimensions
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  this->WriteAsFloat(static_cast<float>(xsize), 10);
  this->WriteAsFloat(static_cast<float>(ysize), 14);
  this->WriteAsFloat(static_cast<float>(zsize), 18);
  this->WriteAsFloat(float(0), 22);

  // Write patient description
  for (i = 0; i < 80; i++) {
    this->WriteAsChar(char(0), 26+i*sizeof(char));
  }

  _input->GetOrientation(xaxis, yaxis, zaxis);
  this->WriteAsFloat(static_cast<float>(xaxis[0]), 106);
  this->WriteAsFloat(static_cast<float>(xaxis[1]), 110);
  this->WriteAsFloat(static_cast<float>(xaxis[2]), 114);
  this->WriteAsFloat(static_cast<float>(yaxis[0]), 122);
  this->WriteAsFloat(static_cast<float>(yaxis[1]), 126);
  this->WriteAsFloat(static_cast<float>(yaxis[2]), 130);

  // Write identifier
  this->WriteAsInt(0, 154);

  // Write spare bytes
  for (i = 0; i < 28; i++) {
    this->WriteAsChar(char(0), 158+i*sizeof(char));
  }

  // Write flag and orientation
  this->WriteAsChar(char(0), 186);
  this->WriteAsChar(char(0), 187);

  // Write min and max values
  _input->GetMinMax(&min, &max);
  this->WriteAsDouble(min, 188);
  this->WriteAsDouble(max, 196);

  // Write origin
  _input->GetOrigin(pos[0], pos[1], pos[2]);
  this->WriteAsDouble(pos[0], 204);
  this->WriteAsDouble(pos[1], 212);
  this->WriteAsDouble(pos[2], 220);
  this->WriteAsDouble(double(0), 228);

  // Write magic number
  this->WriteAsFloat(float(0), 236);
  this->WriteAsFloat(float(0), 240);
  this->WriteAsUInt(GIPL_MAGIC_EXT, 244);
  this->WriteAsFloat(float(0), 248);
  this->WriteAsUInt(GIPL_MAGIC1, 252);

  // Calculate data address
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = GIPL_HEADERSIZE + i * _input->GetX() * _input->GetY() * sizeof(irtkGreyPixel);
  }
}

template <> void irtkImageToFileGIPL<irtkRealPixel>::Initialize()
{
  int i;
  irtkRealPixel min, max;
  double xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3], pos[3];

  // Initialize base class
  this->irtkImageToFile<irtkRealPixel>::Initialize();

  // Write image dimensions
  this->WriteAsShort(_input->GetX(), 0);
  this->WriteAsShort(_input->GetY(), 2);
  this->WriteAsShort(_input->GetZ(), 4);
  this->WriteAsShort(1, 6);

  // Write type
  this->WriteAsShort(GIPL_FLOAT, 8);

  // Write voxel dimensions
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  this->WriteAsFloat(static_cast<float>(xsize), 10);
  this->WriteAsFloat(static_cast<float>(ysize), 14);
  this->WriteAsFloat(static_cast<float>(zsize), 18);
  this->WriteAsFloat(float(0), 22);

  // Write patient description
  for (i = 0; i < 80; i++) {
    this->WriteAsChar(char(0), 26+i*sizeof(char));
  }

  // Write image orientation
  _input->GetOrientation(xaxis, yaxis, zaxis);
  this->WriteAsFloat(static_cast<float>(xaxis[0]), 106);
  this->WriteAsFloat(static_cast<float>(xaxis[1]), 110);
  this->WriteAsFloat(static_cast<float>(xaxis[2]), 114);
  this->WriteAsFloat(static_cast<float>(yaxis[0]), 122);
  this->WriteAsFloat(static_cast<float>(yaxis[1]), 126);
  this->WriteAsFloat(static_cast<float>(yaxis[2]), 130);

  // Write identifier
  this->WriteAsInt(0, 154);

  // Write spare bytes
  for (i = 0; i < 28; i++) {
    this->WriteAsChar(char(0), 158+i*sizeof(char));
  }

  // Write flag and orientation
  this->WriteAsChar(char(0), 186);
  this->WriteAsChar(char(0), 187);

  // Write min and max values
  _input->GetMinMax(&min, &max);
  this->WriteAsDouble(min, 188);
  this->WriteAsDouble(max, 196);

  // Write origin
  _input->GetOrigin(pos[0], pos[1], pos[2]);
  this->WriteAsDouble(pos[0], 204);
  this->WriteAsDouble(pos[1], 212);
  this->WriteAsDouble(pos[2], 220);
  this->WriteAsDouble(double(0), 228);

  // Write magic number
  this->WriteAsFloat(float(0), 236);
  this->WriteAsFloat(float(0), 240);
  this->WriteAsUInt(GIPL_MAGIC_EXT, 244);
  this->WriteAsFloat(float(0), 248);
  this->WriteAsUInt(GIPL_MAGIC1, 252);

  // Calculate data address
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = GIPL_HEADERSIZE + i * _input->GetX() * _input->GetY() * sizeof(irtkRealPixel);
  }
}

template <> void irtkImageToFileGIPL<irtkVector3D<char> >::Initialize()
{
  int i;
  irtkVector3D<char> min, max;
  double xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3], pos[3];

  // Initialize base class
  this->irtkImageToFile<irtkVector3D<char> >::Initialize();

  // Write the iamge dimensions.
  this->WriteAsShort(_input->GetX(), 0);
  this->WriteAsShort(_input->GetY(), 2);
  this->WriteAsShort(_input->GetZ(), 4);
  this->WriteAsShort(1, 6);

  // Write type
  this->WriteAsShort(GIPL_VECTOR_3D_CHAR, 8);

  // Write voxel dimensions
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  this->WriteAsFloat(static_cast<float>(xsize), 10);
  this->WriteAsFloat(static_cast<float>(ysize), 14);
  this->WriteAsFloat(static_cast<float>(zsize), 18);
  this->WriteAsFloat(float(0), 22);

  // Write patient description.
  for (i = 0; i < 80; i++) {
    this->WriteAsChar(char(0), 26+i*sizeof(char));
  }

  // Write image orientation.
  _input->GetOrientation(xaxis, yaxis, zaxis);
  this->WriteAsFloat(static_cast<float>(xaxis[0]), 106);
  this->WriteAsFloat(static_cast<float>(xaxis[1]), 110);
  this->WriteAsFloat(static_cast<float>(xaxis[2]), 114);
  this->WriteAsFloat(static_cast<float>(yaxis[0]), 122);
  this->WriteAsFloat(static_cast<float>(yaxis[1]), 126);
  this->WriteAsFloat(static_cast<float>(yaxis[2]), 130);

  // Write identifier.
  this->WriteAsInt(0, 154);

  // Write spare bytes
  for (i = 0; i < 28; i++) {
    this->WriteAsChar(char(0), 158+i*sizeof(char));
  }

  // Write flag and orientation.
  this->WriteAsChar(char(0), 186);
  this->WriteAsChar(char(0), 187);

  // Write min and max values.
  _input->GetMinMax(&min, &max);
  double lmin, lmax;
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  this->WriteAsDouble(lmin, 188);
  this->WriteAsDouble(lmax, 196);

  // Write origin.
  _input->GetOrigin(pos[0], pos[1], pos[2]);
  this->WriteAsDouble(pos[0], 204);
  this->WriteAsDouble(pos[1], 212);
  this->WriteAsDouble(pos[2], 220);
  this->WriteAsDouble(double(0), 228);

  // Write magic number.
  this->WriteAsFloat(float(0), 236);
  this->WriteAsFloat(float(0), 240);
  this->WriteAsUInt(GIPL_MAGIC_EXT, 244);
  this->WriteAsFloat(float(0), 248);
  this->WriteAsUInt(GIPL_MAGIC1, 252);

  // Calculate data address.
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = GIPL_HEADERSIZE + i*_input->GetX()*_input->GetY()*sizeof(char)*3;
  }
}

template <> void irtkImageToFileGIPL<irtkVector3D<short> >::Initialize()
{
  int i;
  irtkVector3D<short> min, max;
  double xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3], pos[3];

  // Initialize base class
  this->irtkImageToFile<irtkVector3D<short> >::Initialize();

  // Write the iamge dimensions.
  this->WriteAsShort(_input->GetX(), 0);
  this->WriteAsShort(_input->GetY(), 2);
  this->WriteAsShort(_input->GetZ(), 4);
  this->WriteAsShort(1, 6);

  // Write type
  this->WriteAsShort(GIPL_VECTOR_3D_SHORT, 8);

  // Write voxel dimensions
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  this->WriteAsFloat(static_cast<float>(xsize), 10);
  this->WriteAsFloat(static_cast<float>(ysize), 14);
  this->WriteAsFloat(static_cast<float>(zsize), 18);
  this->WriteAsFloat(float(0), 22);

  // Write patient description.
  for (i = 0; i < 80; i++) {
    this->WriteAsChar(char(0), 26+i*sizeof(char));
  }

  // Write image orientation.
  _input->GetOrientation(xaxis, yaxis, zaxis);
  this->WriteAsFloat(static_cast<float>(xaxis[0]), 106);
  this->WriteAsFloat(static_cast<float>(xaxis[1]), 110);
  this->WriteAsFloat(static_cast<float>(xaxis[2]), 114);
  this->WriteAsFloat(static_cast<float>(yaxis[0]), 122);
  this->WriteAsFloat(static_cast<float>(yaxis[1]), 126);
  this->WriteAsFloat(static_cast<float>(yaxis[2]), 130);

  // Write identifier.
  this->WriteAsInt(0, 154);

  // Write spare bytes
  for (i = 0; i < 28; i++) {
    this->WriteAsChar(char(0), 158+i*sizeof(char));
  }

  // Write flag and orientation.
  this->WriteAsChar(char(0), 186);
  this->WriteAsChar(char(0), 187);

  // Write min and max values.
  _input->GetMinMax(&min, &max);
  double lmin, lmax;
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  this->WriteAsDouble(lmin, 188);
  this->WriteAsDouble(lmax, 196);

  // Write origin.
  _input->GetOrigin(pos[0], pos[1], pos[2]);
  this->WriteAsDouble(pos[0], 204);
  this->WriteAsDouble(pos[1], 212);
  this->WriteAsDouble(pos[2], 220);
  this->WriteAsDouble(double(0), 228);

  // Write magic number.
  this->WriteAsFloat(float(0), 236);
  this->WriteAsFloat(float(0), 240);
  this->WriteAsUInt(GIPL_MAGIC_EXT, 244);
  this->WriteAsFloat(float(0), 248);
  this->WriteAsUInt(GIPL_MAGIC1, 252);

  // Calculate data address.
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = GIPL_HEADERSIZE + i*_input->GetX()*_input->GetY()*sizeof(short)*3;
  }
}

template <> void irtkImageToFileGIPL<irtkVector3D<float> >::Initialize()
{
  int i;
  irtkVector3D<float> min, max;
  double xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3], pos[3];

  // Initialize base class
  this->irtkImageToFile<irtkVector3D<float> >::Initialize();

  // Write the iamge dimensions.
  this->WriteAsShort(_input->GetX(), 0);
  this->WriteAsShort(_input->GetY(), 2);
  this->WriteAsShort(_input->GetZ(), 4);
  this->WriteAsShort(1, 6);

  // Write type
  this->WriteAsShort(GIPL_VECTOR_3D_FLOAT, 8);

  // Write voxel dimensions
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  this->WriteAsFloat(static_cast<float>(xsize), 10);
  this->WriteAsFloat(static_cast<float>(ysize), 14);
  this->WriteAsFloat(static_cast<float>(zsize), 18);
  this->WriteAsFloat(float(0), 22);

  // Write patient description.
  for (i = 0; i < 80; i++) {
    this->WriteAsChar(char(0), 26+i*sizeof(char));
  }

  // Write image orientation.
  _input->GetOrientation(xaxis, yaxis, zaxis);
  this->WriteAsFloat(static_cast<float>(xaxis[0]), 106);
  this->WriteAsFloat(static_cast<float>(xaxis[1]), 110);
  this->WriteAsFloat(static_cast<float>(xaxis[2]), 114);
  this->WriteAsFloat(static_cast<float>(yaxis[0]), 122);
  this->WriteAsFloat(static_cast<float>(yaxis[1]), 126);
  this->WriteAsFloat(static_cast<float>(yaxis[2]), 130);

  // Write identifier.
  this->WriteAsInt(0, 154);

  // Write spare bytes
  for (i = 0; i < 28; i++) {
    this->WriteAsChar(char(0), 158+i*sizeof(char));
  }

  // Write flag and orientation.
  this->WriteAsChar(char(0), 186);
  this->WriteAsChar(char(0), 187);

  // Write min and max values.
  _input->GetMinMax(&min, &max);
  double lmin, lmax;
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  this->WriteAsDouble(lmin, 188);
  this->WriteAsDouble(lmax, 196);

  // Write origin.
  _input->GetOrigin(pos[0], pos[1], pos[2]);
  this->WriteAsDouble(pos[0], 204);
  this->WriteAsDouble(pos[1], 212);
  this->WriteAsDouble(pos[2], 220);
  this->WriteAsDouble(double(0), 228);

  // Write magic number.
  this->WriteAsFloat(float(0), 236);
  this->WriteAsFloat(float(0), 240);
  this->WriteAsUInt(GIPL_MAGIC_EXT, 244);
  this->WriteAsFloat(float(0), 248);
  this->WriteAsUInt(GIPL_MAGIC1, 252);

  // Calculate data address.
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = GIPL_HEADERSIZE + i*_input->GetX()*_input->GetY()*sizeof(float)*3;
  }
}

template <> void irtkImageToFileGIPL<irtkVector3D<double> >::Initialize()
{
  int i;
  irtkVector3D<double> min, max;
  double xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3], pos[3];

  // Initialize base class
  this->irtkImageToFile<irtkVector3D<double> >::Initialize();

  // Write the iamge dimensions.
  this->WriteAsShort(_input->GetX(), 0);
  this->WriteAsShort(_input->GetY(), 2);
  this->WriteAsShort(_input->GetZ(), 4);
  this->WriteAsShort(1, 6);

  // Write type
  this->WriteAsShort(GIPL_VECTOR_3D_DOUBLE, 8);

  // Write voxel dimensions
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  this->WriteAsFloat(static_cast<float>(xsize), 10);
  this->WriteAsFloat(static_cast<float>(ysize), 14);
  this->WriteAsFloat(static_cast<float>(zsize), 18);
  this->WriteAsFloat(float(0), 22);

  // Write patient description.
  for (i = 0; i < 80; i++) {
    this->WriteAsChar(char(0), 26+i*sizeof(char));
  }

  // Write image orientation.
  _input->GetOrientation(xaxis, yaxis, zaxis);
  this->WriteAsFloat(static_cast<float>(xaxis[0]), 106);
  this->WriteAsFloat(static_cast<float>(xaxis[1]), 110);
  this->WriteAsFloat(static_cast<float>(xaxis[2]), 114);
  this->WriteAsFloat(static_cast<float>(yaxis[0]), 122);
  this->WriteAsFloat(static_cast<float>(yaxis[1]), 126);
  this->WriteAsFloat(static_cast<float>(yaxis[2]), 130);

  // Write identifier.
  this->WriteAsInt(0, 154);

  // Write spare bytes
  for (i = 0; i < 28; i++) {
    this->WriteAsChar(char(0), 158+i*sizeof(char));
  }

  // Write flag and orientation.
  this->WriteAsChar(char(0), 186);
  this->WriteAsChar(char(0), 187);

  // Write min and max values.
  _input->GetMinMax(&min, &max);
  double lmin, lmax;
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  this->WriteAsDouble(lmin, 188);
  this->WriteAsDouble(lmax, 196);

  // Write origin.
  _input->GetOrigin(pos[0], pos[1], pos[2]);
  this->WriteAsDouble(pos[0], 204);
  this->WriteAsDouble(pos[1], 212);
  this->WriteAsDouble(pos[2], 220);
  this->WriteAsDouble(double(0), 228);

  // Write magic number.
  this->WriteAsFloat(float(0), 236);
  this->WriteAsFloat(float(0), 240);
  this->WriteAsUInt(GIPL_MAGIC_EXT, 244);
  this->WriteAsFloat(float(0), 248);
  this->WriteAsUInt(GIPL_MAGIC1, 252);

  // Calculate data address.
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = GIPL_HEADERSIZE + i*_input->GetX()*_input->GetY()*sizeof(double)*3;
  }
}

template <> void irtkImageToFileGIPL<irtkRGBPixel>::Initialize()
{
  cerr << "irtkImageToFileGIPL<irtkRGBPixel>::Initialize: Not supported" << endl;
  exit(1);
}

template class irtkImageToFileGIPL<irtkBytePixel>;
template class irtkImageToFileGIPL<irtkGreyPixel>;
template class irtkImageToFileGIPL<irtkRealPixel>;
template class irtkImageToFileGIPL<irtkRGBPixel>;
template class irtkImageToFileGIPL<irtkVector3D<char> >;
template class irtkImageToFileGIPL<irtkVector3D<short> >;
template class irtkImageToFileGIPL<irtkVector3D<float> >;
template class irtkImageToFileGIPL<irtkVector3D<double> >;


