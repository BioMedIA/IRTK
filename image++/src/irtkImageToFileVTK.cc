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

#include <irtkVTK.h>

template <class VoxelType> const char *irtkImageToFileVTK<VoxelType>::NameOfClass()
{
  return "irtkImageToFileVTK";
}

template <> void irtkImageToFileVTK<irtkBytePixel>::Initialize()
{
  irtkPoint p;
  int i, x, y, z;
  char header[2048];
  double xsize, ysize, zsize;

  // Initialize base class
  this->irtkImageToFile<irtkBytePixel>::Initialize();

  // Get image information from input
  x = _input->GetX();
  y = _input->GetY();
  z = _input->GetZ();

  double xaxis[3] = { 1, 0, 0 };
  double yaxis[3] = { 0, 1, 0 };
  double zaxis[3] = { 0, 0, 1 };

  _input->PutOrientation(xaxis, yaxis, zaxis);

  _input->ImageToWorld(p);
  _input->GetPixelSize(&xsize, &ysize, &zsize);

  // Construct header
  sprintf(header, "%s\nCreated by %s\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\nSCALARS scalars unsigned_char 1\nLOOKUP_TABLE default\n", VTK_MAGIC2, this->NameOfClass(), x, y, z, p._x, p._y, p._z, xsize, ysize, zsize, x*y*z);

  // Write magic number
  this->Write(header, 0, strlen(header));

  // Calculate data address
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = strlen(header) + i * _input->GetX() * _input->GetY() * sizeof(irtkBytePixel);
  }
}

template <> void irtkImageToFileVTK<irtkGreyPixel>::Initialize()
{
  irtkPoint p;
  int i, x, y, z;
  char header[2048];
  double xsize, ysize, zsize;

  // Initialize base class
  this->irtkImageToFile<irtkGreyPixel>::Initialize();

  // Get image information from input
  x = _input->GetX();
  y = _input->GetY();
  z = _input->GetZ();

  double xaxis[3] = { 1, 0, 0 };
  double yaxis[3] = { 0, 1, 0 };
  double zaxis[3] = { 0, 0, 1 };

  _input->PutOrientation(xaxis, yaxis, zaxis);

  _input->ImageToWorld(p);
  _input->GetPixelSize(&xsize, &ysize, &zsize);
  // Construct header
  sprintf(header, "%s\nCreated by %s\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\nSCALARS scalars short 1\nLOOKUP_TABLE default\n", VTK_MAGIC2, this->NameOfClass(), x, y, z, p._x, p._y, p._z, xsize, ysize, zsize, x*y*z);

  // Write magic number
  this->Write(header, 0, strlen(header));

  // Calculate data address
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = strlen(header) + i * _input->GetX() * _input->GetY() * sizeof(irtkGreyPixel);
  }
}

template <> void irtkImageToFileVTK<irtkRealPixel>::Initialize()
{
  irtkPoint p;
  int i, x, y, z;
  char header[2048];
  double xsize, ysize, zsize;

  // Initialize base class
  this->irtkImageToFile<irtkRealPixel>::Initialize();

  // Get image information from input
  x = _input->GetX();
  y = _input->GetY();
  z = _input->GetZ();

  double xaxis[3] = { 1, 0, 0 };
  double yaxis[3] = { 0, 1, 0 };
  double zaxis[3] = { 0, 0, 1 };

  _input->PutOrientation(xaxis, yaxis, zaxis);

  _input->ImageToWorld(p);
  _input->GetPixelSize(&xsize, &ysize, &zsize);

  // Construct header
  sprintf(header, "%s\nCreated by %s\nBINARY\nDATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\nORIGIN %f %f %f\nSPACING %f %f %f\nPOINT_DATA %d\nSCALARS scalars float 1\nLOOKUP_TABLE default\n", VTK_MAGIC2, this->NameOfClass(), x, y, z, p._x, p._y, p._z, xsize, ysize, zsize, x*y*z);

  // Write magic number
  this->Write(header, 0, strlen(header));

  // Calculate data address
  for (i = 0; i < _input->GetZ(); i++) {
    _addr[i] = strlen(header) + i * _input->GetX() * _input->GetY() * sizeof(irtkRealPixel);
  }
}

template <> void irtkImageToFileVTK<irtkRGBPixel>::Initialize()
{
  cerr << "irtkImageToFileVTK<irtkRGBPixel>::Initialize: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileVTK<irtkVector3D<char> >::Initialize()
{
  cerr << "irtkImageToFileVTK<irtkVector3D<char> >::Initialize: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileVTK<irtkVector3D<short> >::Initialize()
{
  cerr << "irtkImageToFileVTK<irtkVector3D<short> >::Initialize: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileVTK<irtkVector3D<float> >::Initialize()
{
  cerr << "irtkImageToFileVTK<irtkVector3D<float> >::Initialize: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileVTK<irtkVector3D<double> >::Initialize()
{
  cerr << "irtkImageToFileVTK<irtkVector3D<double> >::Initialize: Not supported" << endl;
  exit(1);
}

template class irtkImageToFileVTK<irtkBytePixel>;
template class irtkImageToFileVTK<irtkGreyPixel>;
template class irtkImageToFileVTK<irtkRealPixel>;
template class irtkImageToFileVTK<irtkRGBPixel>;
template class irtkImageToFileVTK<irtkVector3D<char> >;
template class irtkImageToFileVTK<irtkVector3D<short> >;
template class irtkImageToFileVTK<irtkVector3D<float> >;
template class irtkImageToFileVTK<irtkVector3D<double> >;
