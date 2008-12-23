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

#include <irtkANALYZE.h>

template <class VoxelType> irtkImageToFileANALYZE<VoxelType>::irtkImageToFileANALYZE() : irtkImageToFile<VoxelType>()
{
  this->_headername = NULL;
  this->_reflectY = True;
}

template <class VoxelType> irtkImageToFileANALYZE<VoxelType>::~irtkImageToFileANALYZE()
{
  if (this->_headername != NULL) free(this->_headername);
}

template <class VoxelType> void irtkImageToFileANALYZE<VoxelType>::SetOutput(const char *name)
{
  int length;
  char *imagename;

  if (this->_headername != NULL) free(this->_headername);
  if (strstr(name, ".gz") == NULL) {
    this->_headername = strdup(name);
    imagename   = strdup(name);
    length      = strlen(name);
    imagename[length-1] = 'g';
    imagename[length-2] = 'm';
    imagename[length-3] = 'i';
  } else {
    this->_headername = strdup(name);
    imagename   = strdup(name);
    length      = strlen(name);
    imagename[length-4] = 'g';
    imagename[length-5] = 'm';
    imagename[length-6] = 'i';
  }
  this->irtkImageToFile<VoxelType>::SetOutput(imagename);
}

template <class VoxelType> const char *irtkImageToFileANALYZE<VoxelType>::NameOfClass()
{
  return "irtkImageToFileANALYZE";
}

template <> void irtkImageToFileANALYZE<irtkBytePixel>::Initialize()
{
  int i, x, y, z;
  irtkCofstream to;
  irtkBytePixel min, max;
  double xsize, ysize, zsize;
  irtkANALYZEHeader hdr;

  x = this->_input->GetX();
  y = this->_input->GetY();
  z = this->_input->GetZ();
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);
  this->_input->GetMinMax(&min, &max);

  hdr.bits       = 8;
  hdr.data_type  = ANALYZE_UNSIGNED_CHAR;
  hdr.dims[0]    = 4;
  hdr.dims[1]    = this->_input->GetX();
  hdr.dims[2]    = this->_input->GetY();
  hdr.dims[3]    = this->_input->GetZ();
  hdr.pixdims[0] = 1;
  hdr.pixdims[1] = xsize;
  hdr.pixdims[2] = ysize;
  hdr.pixdims[3] = zsize;
  hdr.glmax      = max;
  hdr.glmin      = min;

  // Write header
  hdr.Write(this->_headername);

  // Initialize base class
  this->irtkImageToFile<irtkBytePixel>::Initialize();

  // Calculate data address
  for (i = 0; i < this->_input->GetZ(); i++) {
    this->_addr[i] = i * this->_input->GetX() * this->_input->GetY() * sizeof(irtkBytePixel);
  }
}

template <> void irtkImageToFileANALYZE<irtkGreyPixel>::Initialize()
{
  int i, x, y, z;
  irtkCofstream to;
  irtkGreyPixel min, max;
  double xsize, ysize, zsize;
  irtkANALYZEHeader hdr;

  x = this->_input->GetX();
  y = this->_input->GetY();
  z = this->_input->GetZ();
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);
  this->_input->GetMinMax(&min, &max);

  hdr.bits       = 16;
  hdr.data_type  = ANALYZE_SIGNED_SHORT;
  hdr.dims[0]    = 4;
  hdr.dims[1]    = this->_input->GetX();
  hdr.dims[2]    = this->_input->GetY();
  hdr.dims[3]    = this->_input->GetZ();
  hdr.pixdims[0] = 1;
  hdr.pixdims[1] = xsize;
  hdr.pixdims[2] = ysize;
  hdr.pixdims[3] = zsize;
  hdr.glmax      = max;
  hdr.glmin      = min;

  // Write header
  hdr.Write(this->_headername);

  // Initialize base class
  this->irtkImageToFile<irtkGreyPixel>::Initialize();

  // Calculate data address
  for (i = 0; i < this->_input->GetZ(); i++) {
    this->_addr[i] = i * this->_input->GetX() * this->_input->GetY() * sizeof(irtkGreyPixel);
  }
}

template <> void irtkImageToFileANALYZE<irtkRealPixel>::Initialize()
{
  int i, x, y, z;
  irtkCofstream to;
  irtkRealPixel min, max;
  double xsize, ysize, zsize;
  irtkANALYZEHeader hdr;

  x = this->_input->GetX();
  y = this->_input->GetY();
  z = this->_input->GetZ();
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);
  this->_input->GetMinMax(&min, &max);

  hdr.bits       = 32;
  hdr.data_type  = ANALYZE_FLOAT;
  hdr.dims[0]    = 4;
  hdr.dims[1]    = this->_input->GetX();
  hdr.dims[2]    = this->_input->GetY();
  hdr.dims[3]    = this->_input->GetZ();
  hdr.pixdims[0] = 1;
  hdr.pixdims[1] = xsize;
  hdr.pixdims[2] = ysize;
  hdr.pixdims[3] = zsize;
  hdr.glmax      = round(max);
  hdr.glmin      = round(min);

  // Write header
  hdr.Write(this->_headername);

  // Initialize base class
  this->irtkImageToFile<irtkRealPixel>::Initialize();

  // Calculate data address
  for (i = 0; i < this->_input->GetZ(); i++) {
    this->_addr[i] = i * this->_input->GetX() * this->_input->GetY() * sizeof(irtkRealPixel);
  }
}

template <> void irtkImageToFileANALYZE<irtkRGBPixel>::Initialize()
{
  cerr << "irtkImageToFileANALYZE<irtkRGBPixel>::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileANALYZE<irtkVector3D<char> >::Initialize()
{
  cerr << "irtkImageToFileANALYZE<irtkVector3D<char> >::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileANALYZE<irtkVector3D<short> >::Initialize()
{
  cerr << "irtkImageToFileANALYZE<irtkVector3D<short> >::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileANALYZE<irtkVector3D<float> >::Initialize()
{
  cerr << "irtkImageToFileANALYZE<irtkVector3D<float> >::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileANALYZE<irtkVector3D<double> >::Initialize()
{
  cerr << "irtkImageToFileANALYZE<irtkVector3D<double> >::Run: Not supported" << endl;
  exit(1);
}

template class irtkImageToFileANALYZE<irtkBytePixel>;
template class irtkImageToFileANALYZE<irtkGreyPixel>;
template class irtkImageToFileANALYZE<irtkRealPixel>;
template class irtkImageToFileANALYZE<irtkRGBPixel>;
template class irtkImageToFileANALYZE<irtkVector3D<char> >;
template class irtkImageToFileANALYZE<irtkVector3D<short> >;
template class irtkImageToFileANALYZE<irtkVector3D<float> >;
template class irtkImageToFileANALYZE<irtkVector3D<double> >;

