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

irtkImageToFileANALYZE::irtkImageToFileANALYZE() : irtkImageToFile()
{
  this->_headername = NULL;
  this->_reflectY = true;
}

irtkImageToFileANALYZE::~irtkImageToFileANALYZE()
{
  if (this->_headername != NULL) free(this->_headername);
}

void irtkImageToFileANALYZE::SetOutput(const char *name)
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
  this->irtkImageToFile::SetOutput(imagename);
}

const char *irtkImageToFileANALYZE::NameOfClass()
{
  return "irtkImageToFileANALYZE";
}

void irtkImageToFileANALYZE::Initialize()
{
  int x, y, z;
  irtkCofstream to;
  double min, max;
  double xsize, ysize, zsize;
  irtkANALYZEHeader hdr;

  x = this->_input->GetX();
  y = this->_input->GetY();
  z = this->_input->GetZ();
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);
  this->_input->GetMinMaxAsDouble(&min, &max);

  switch (this->_input->GetScalarType()) {
    case IRTK_VOXEL_UNSIGNED_CHAR: {
        hdr.bits       = 8;
        hdr.data_type  = ANALYZE_UNSIGNED_CHAR;
        break;
      }
    case IRTK_VOXEL_SHORT: {
        hdr.bits       = 16;
        hdr.data_type  = ANALYZE_SIGNED_SHORT;
        break;
      }
   case IRTK_VOXEL_FLOAT: {
        hdr.bits       = 32;
        hdr.data_type  = ANALYZE_FLOAT;
        break;
      }
    case IRTK_VOXEL_DOUBLE: {
        hdr.bits       = 64;
        hdr.data_type  = ANALYZE_DOUBLE;
        break;
      }
    default:
      cerr << "irtkImageToFileANALYZE::Initialize: Not supported for this image type" << endl;
      exit(1);
  }

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
  this->irtkImageToFile::Initialize();

  // Calculate data address
  _start = 0;
}

