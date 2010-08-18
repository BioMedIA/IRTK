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

#include <irtkGIPL.h>

const char *irtkFileGIPLToImage::NameOfClass()
{
  return "irtkFileGIPLToImage";
}

int irtkFileGIPLToImage::CheckHeader(const char *filename)
{
  unsigned int magic_number;

  // Create file stream
  irtkCifstream from;

  // Open new file for reading
  from.Open(filename);

  // Read header
  from.ReadAsUInt(&magic_number, 1, 252);

  // Close file
  from.Close();

  // Check header
  if ((magic_number != GIPL_MAGIC1) && (magic_number != GIPL_MAGIC2)) {
    return False;
  } else {
    return True;
  }
}

void irtkFileGIPLToImage::ReadHeader()
{
  float size;
  double origin;
  unsigned short dim, type;
  unsigned int magic_number, magic_ext;

  // Read header
  this->ReadAsUInt(&magic_number, 1, 252);

  // Check header
  if ((magic_number != GIPL_MAGIC1) && (magic_number != GIPL_MAGIC2)) {
    cerr << this->NameOfClass() << "::Read_Header: Can't read magic number"
         << endl;
    exit(1);
  }

  // Read image dimensions
  this->ReadAsUShort(&dim, 1, 0);
  this->_attr._x = dim;
  this->ReadAsUShort(&dim, 1, 2);
  this->_attr._y = dim;
  this->ReadAsUShort(&dim, 1, 4);
  this->_attr._z = dim;
  this->ReadAsUShort(&dim, 1, 6);
  this->_attr._t = dim;
  if(this->_attr._t <= 0)
	  this->_attr._t = 1;

  // Read voxel dimensions
  this->ReadAsFloat(&size, 1, 10);
  this->_attr._dx = size;
  this->ReadAsFloat(&size, 1, 14);
  this->_attr._dy = size;
  this->ReadAsFloat(&size, 1, 18);
  this->_attr._dz = size;
  this->ReadAsFloat(&size, 1, 22);
  this->_attr._dt = size;
  if(this->_attr._dt <= 0)
	  this->_attr._dt = 1;

  // Read extension flag
  this->ReadAsUInt(&magic_ext, 1, 244);

  // Check if additional information about orientation and origin is available
  if (magic_ext == GIPL_MAGIC_EXT) {

    // Read image orientation if available
    this->ReadAsFloat(&size, 1, 106);
    this->_attr._xaxis[0] = size;
    this->ReadAsFloat(&size, 1, 110);
    this->_attr._xaxis[1] = size;
    this->ReadAsFloat(&size, 1, 114);
    this->_attr._xaxis[2] = size;
    this->ReadAsFloat(&size, 1, 122);
    this->_attr._yaxis[0] = size;
    this->ReadAsFloat(&size, 1, 126);
    this->_attr._yaxis[1] = size;
    this->ReadAsFloat(&size, 1, 130);
    this->_attr._yaxis[2] = size;

    // Construct the z-axis.
    // Assume a right handed (`neurological') set of axes!
    this->_attr._zaxis[0] = this->_attr._xaxis[1]*this->_attr._yaxis[2] - this->_attr._xaxis[2]*this->_attr._yaxis[1];
    this->_attr._zaxis[1] = this->_attr._xaxis[2]*this->_attr._yaxis[0] - this->_attr._xaxis[0]*this->_attr._yaxis[2];
    this->_attr._zaxis[2] = this->_attr._xaxis[0]*this->_attr._yaxis[1] - this->_attr._xaxis[1]*this->_attr._yaxis[0];

    // Read image origin if available
    this->ReadAsDouble(&origin, 1, 204);
    this->_attr._xorigin = origin;
    this->ReadAsDouble(&origin, 1, 212);
    this->_attr._yorigin = origin;
    this->ReadAsDouble(&origin, 1, 220);
    this->_attr._zorigin = origin;
  }

  // Read type of voxels
  this->ReadAsUShort(&type, 1, 8);

  // Calculate type of voxels and number of bytes per voxel
  switch (type) {
  case GIPL_CHAR:
    this->_type  = IRTK_VOXEL_CHAR;
    this->_bytes = 1;
    break;
  case GIPL_U_CHAR:
    this->_type  = IRTK_VOXEL_UNSIGNED_CHAR;
    this->_bytes = 1;
    break;
  case GIPL_SHORT:
    this->_type  = IRTK_VOXEL_SHORT;
    this->_bytes = 2;
    break;
  case GIPL_U_SHORT:
    this->_type  = IRTK_VOXEL_UNSIGNED_SHORT;
    this->_bytes = 2;
    break;
  case GIPL_INT:
    this->_type  = IRTK_VOXEL_INT;
    this->_bytes = 4;
    break;
  case GIPL_U_INT:
    this->_type  = IRTK_VOXEL_UNSIGNED_INT;
    this->_bytes = 4;
    break;
  case GIPL_FLOAT:
    this->_type  = IRTK_VOXEL_FLOAT;
    this->_bytes = 4;
    break;
  case GIPL_DOUBLE:
    this->_type  = IRTK_VOXEL_DOUBLE;
    this->_bytes = 8;
    break;
  default:
    cerr << this->NameOfClass() << "::Read_Header(): Unknown voxel type"
         << endl;
    exit(1);
  }

  // Data starts here
  this->_start = GIPL_HEADERSIZE;
}