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

#include <irtkFileToImage.h>

#include <irtkANALYZE.h>

#include <sys/types.h>
#include <sys/stat.h>

irtkFileANALYZEToImage::irtkFileANALYZEToImage()
{
  _headername = NULL;

  // Analyze-specific
  this->_reflectY = true;
}

irtkFileANALYZEToImage::~irtkFileANALYZEToImage()
{
  if (this->_headername != NULL) free(this->_headername);
}

const char *irtkFileANALYZEToImage::NameOfClass()
{
  return "irtkFileANALYZEToImage";
}

int irtkFileANALYZEToImage::CheckHeader(const char *filename)
{
  if ((strstr(basename2(filename), ".hdr") == NULL) && (strstr(basename2(filename), ".HDR") == NULL)) {
    return false;
  } else {
    return true;
  }
}

void irtkFileANALYZEToImage::SetInput(const char *headername)
{
  int length;
  char imagename[255];
  struct stat buf;

  // Delete old file name
  if (this->_headername != NULL) free(this->_headername);

  // Copy new file name
  this->_headername = strdup(headername);

  // Parse img file name
  if (strstr(basename2(headername), ".hdr.gz") != NULL) {
    length = strlen(headername);
    sprintf(imagename, "%s", headername);
    imagename[length-1] = 'z';
    imagename[length-2] = 'g';
    imagename[length-3] = '.';
    imagename[length-4] = 'g';
    imagename[length-5] = 'm';
    imagename[length-6] = 'i';

    // Check if compressed file exists
    if (stat(imagename, &buf) != 0) {
      cerr << this->NameOfClass() << ": Can't open file " << imagename << endl;
      exit(1);
    }
  } else {
    length = strlen(headername);
    sprintf(imagename, "%s", headername);
    imagename[length-1] = 'g';
    imagename[length-2] = 'm';
    imagename[length-3] = 'i';

    // Check if uncompressed file exists
    if (stat(imagename, &buf) != 0) {
      sprintf(imagename, "%s.gz", headername);
      imagename[length-1] = 'g';
      imagename[length-2] = 'm';
      imagename[length-3] = 'i';

      // Check if gzip compressed file exists
      if (stat(imagename, &buf) != 0) {
        sprintf(imagename, "%s.Z", headername);
        imagename[length-1] = 'g';
        imagename[length-2] = 'm';
        imagename[length-3] = 'i';
        if (stat(imagename, &buf) != 0) {
          cerr << this->NameOfClass() << ": Can't open file " << imagename << endl;
          exit(1);
        }
      }
    }
  }

  this->irtkFileToImage::SetInput(imagename);
}

void irtkFileANALYZEToImage::ReadHeader()
{
  irtkANALYZEHeader hdr;

  // Read header
  hdr.Read(this->_headername);

  // Copy header information
  if (hdr.dims[0] == 3) {
    this->_attr._x  = hdr.dims[1];
    this->_attr._y  = hdr.dims[2];
    this->_attr._z  = hdr.dims[3];
    this->_attr._t  = 1;
    this->_attr._dx = hdr.pixdims[1];
    this->_attr._dy = hdr.pixdims[2];
    this->_attr._dy = hdr.pixdims[3];
    this->_attr._dt = 1;
  } else {
    if (hdr.dims[0] == 4) {
      this->_attr._x  = hdr.dims[1];
      this->_attr._y  = hdr.dims[2];
      this->_attr._z  = hdr.dims[3];
      this->_attr._dx = hdr.pixdims[1];
      this->_attr._dy = hdr.pixdims[2];
      this->_attr._dz = hdr.pixdims[3];
      if (hdr.dims[4] <= 0) {
        this->_attr._t  = 1;
        this->_attr._dt = 1;
      } else {
        this->_attr._t  = hdr.dims[4];
        if (hdr.pixdims[4] < 0.0001) hdr.pixdims[4] = 1;
        this->_attr._dt = hdr.pixdims[4];
      }
    } else {
      cerr << "irtkFileANALYZEToImage::ReadHeader: Invalid no. of dimensions in image" << endl;
      exit(1);
    }
  }

  // If the size of header is not 348 we need to invert swapping
  if (hdr.sizeof_hdr != 348) {
    this->_swapped = !this->_swapped;
  }

  switch (hdr.data_type) {
  case ANALYZE_UNSIGNED_CHAR:
    this->_type  = IRTK_VOXEL_UNSIGNED_CHAR;
    this->_bytes = 1;
    break;
  case ANALYZE_SIGNED_SHORT:
    this->_type  = IRTK_VOXEL_SHORT;
    this->_bytes = 2;
    break;
  case ANALYZE_FLOAT:
    this->_type  = IRTK_VOXEL_FLOAT;
    this->_bytes = 4;
    break;
  case ANALYZE_DOUBLE:
    this->_type  = IRTK_VOXEL_DOUBLE;
    this->_bytes = 8;
    break;
  default:
    cerr << this->NameOfClass() << ": Data type " << hdr.data_type << " not supported, trying signed short data type" << endl;
    this->_type  = IRTK_VOXEL_SHORT;
    this->_bytes = 2;
  }

  // Data starts at 0
  this->_start = 0;
}
