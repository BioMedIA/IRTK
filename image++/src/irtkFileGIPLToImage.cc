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

template <class VoxelType> const char *irtkFileGIPLToImage<VoxelType>::NameOfClass()
{
  return "irtkFileGIPLToImage";
}

template <class VoxelType> int irtkFileGIPLToImage<VoxelType>::CheckHeader(const char *filename)
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

template <class VoxelType> void irtkFileGIPLToImage<VoxelType>::ReadHeader()
{
  int i;
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
  this->_x = dim;
  this->ReadAsUShort(&dim, 1, 2);
  this->_y = dim;
  this->ReadAsUShort(&dim, 1, 4);
  this->_z = dim;
  this->_t = 1;

  // Read voxel dimensions
  this->ReadAsFloat(&size, 1, 10);
  this->_xsize = size;
  this->ReadAsFloat(&size, 1, 14);
  this->_ysize = size;
  this->ReadAsFloat(&size, 1, 18);
  this->_zsize = size;
  this->_tsize = 1;

  // Read extension flag
  this->ReadAsUInt(&magic_ext, 1, 244);

  // Check if additional information about orientation and origin is available
  if (magic_ext == GIPL_MAGIC_EXT) {

    // Read image orientation if available
    this->ReadAsFloat(&size, 1, 106);
    this->_xaxis[0] = size;
    this->ReadAsFloat(&size, 1, 110);
    this->_xaxis[1] = size;
    this->ReadAsFloat(&size, 1, 114);
    this->_xaxis[2] = size;
    this->ReadAsFloat(&size, 1, 122);
    this->_yaxis[0] = size;
    this->ReadAsFloat(&size, 1, 126);
    this->_yaxis[1] = size;
    this->ReadAsFloat(&size, 1, 130);
    this->_yaxis[2] = size;

    // Construct the z-axis.
    // Assume a right handed (`neurological') set of axes!
    this->_zaxis[0] = this->_xaxis[1]*this->_yaxis[2] - this->_xaxis[2]*this->_yaxis[1];
    this->_zaxis[1] = this->_xaxis[2]*this->_yaxis[0] - this->_xaxis[0]*this->_yaxis[2];
    this->_zaxis[2] = this->_xaxis[0]*this->_yaxis[1] - this->_xaxis[1]*this->_yaxis[0];

    // Read image origin if available
    this->ReadAsDouble(&origin, 1, 204);
    this->_xorigin = origin;
    this->ReadAsDouble(&origin, 1, 212);
    this->_yorigin = origin;
    this->ReadAsDouble(&origin, 1, 220);
    this->_zorigin = origin;
  }

  // Read type of voxels
  this->ReadAsUShort(&type, 1, 8);

  // Calculate type of voxels and number of bytes per voxel
  switch (type) {
  case GIPL_CHAR:
    this->_type  = VOXEL_CHAR;
    this->_bytes = 1;
    break;
  case GIPL_U_CHAR:
    this->_type  = VOXEL_U_CHAR;
    this->_bytes = 1;
    break;
  case GIPL_SHORT:
    this->_type  = VOXEL_SHORT;
    this->_bytes = 2;
    break;
  case GIPL_U_SHORT:
    this->_type  = VOXEL_U_SHORT;
    this->_bytes = 2;
    break;
  case GIPL_INT:
    this->_type  = VOXEL_INT;
    this->_bytes = 4;
    break;
  case GIPL_U_INT:
    this->_type  = VOXEL_U_INT;
    this->_bytes = 4;
    break;
  case GIPL_FLOAT:
    this->_type  = VOXEL_FLOAT;
    this->_bytes = 4;
    break;
  case GIPL_VECTOR_3D_CHAR:
    this->_type = VOXEL_VECTOR_3D_CHAR;
    this->_bytes = 3;
    break;
  case GIPL_VECTOR_3D_SHORT:
    this->_type = VOXEL_VECTOR_3D_SHORT;
    this->_bytes = 6;
    break;
  case GIPL_VECTOR_3D_FLOAT:
    this->_type = VOXEL_VECTOR_3D_FLOAT;
    this->_bytes = 12;
    break;
  case GIPL_VECTOR_3D_DOUBLE:
    this->_type = VOXEL_VECTOR_3D_DOUBLE;
    this->_bytes = 24;
    break;
  default:
    cerr << this->NameOfClass() << "::Read_Header(): Unknown voxel type"
         << endl;
    exit(1);
  }

  // Allocate memory for data address
  if (this->_addr != 0) delete []this->_addr;
  this->_addr = new int[this->_z];

  // Calculate data address
  for (i = 0; i < this->_z; i++) {
    this->_addr[i] = GIPL_HEADERSIZE + i * this->_x * this->_y * this->_bytes;
  }
}

template class irtkFileGIPLToImage<irtkBytePixel>;
template class irtkFileGIPLToImage<irtkGreyPixel>;
template class irtkFileGIPLToImage<irtkRealPixel>;
template class irtkFileGIPLToImage<irtkVector3D<char> >;
template class irtkFileGIPLToImage<irtkVector3D<short> >;
template class irtkFileGIPLToImage<irtkVector3D<float> >;
template class irtkFileGIPLToImage<irtkVector3D<double> >;

