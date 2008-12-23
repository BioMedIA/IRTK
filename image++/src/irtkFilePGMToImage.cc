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

template <class VoxelType> const char *irtkFilePGMToImage<VoxelType>::NameOfClass()
{
  return "irtkFilePGMToImage";
}

template <class VoxelType> int irtkFilePGMToImage<VoxelType>::CheckHeader(const char *filename)
{
  char buffer[255];

  // Create file stream
  ifstream from;

  // Open new file for reading
  from.open(filename);

  // Read header, skip comments
  do {
    from.get(buffer, 255);
    from.seekg(1, ios::cur);
  } while (buffer[0] == '#');

  // Close file
  from.close();

  // Check header
  if (strcmp(buffer, PGM_MAGIC) != 0) {
    return False;
  } else {
    return True;
  }
}

template <class VoxelType> void irtkFilePGMToImage<VoxelType>::ReadHeader()
{
  int i;
  char buffer[255];

  // Read header, skip comments
  do {
    this->ReadAsString(buffer, 255);
  } while (buffer[0] == '#');

  // Check header
  if (strcmp(buffer, PGM_MAGIC) != 0) {
    cerr << this->NameOfClass() << "::Read_Header: Can't read magic number: "
         << buffer << endl;
    exit(1);
  }

  // Read voxel dimensions, skip comments
  do {
    this->ReadAsString(buffer, 255);
  } while (buffer[0] == '#');

  // Parse voxel dimensions
  sscanf(buffer, "%d %d", &this->_x, &this->_y);

  // Ignore maximum greyvalue, skip comments
  do {
    this->ReadAsString(buffer, 255);
  } while (buffer[0] == '#');

  // PGM files support only 2D images, so set z and t to 1
  this->_z = 1;
  this->_t = 1;

  // PGM files do not have voxel dimensions, so set them to default values
  this->_xsize = 1;
  this->_ysize = 1;
  this->_zsize = 1;
  this->_tsize = 1;

  // PGM files have voxels which are unsigned char
  this->_type  = VOXEL_U_CHAR;
  this->_bytes = 1;

  // Allocate memory for data address
  if (this->_addr != 0) delete []this->_addr;
  this->_addr = new int[this->_z];

  // Calculate data address
  for (i = 0; i < this->_z; i++) {
    this->_addr[i] = int(this->Tell()) + i * this->_x * this->_y * this->_bytes;
  }
}

template class irtkFilePGMToImage<irtkBytePixel>;
template class irtkFilePGMToImage<irtkGreyPixel>;
template class irtkFilePGMToImage<irtkRealPixel>;
template class irtkFilePGMToImage<irtkVector3D<char> >;
template class irtkFilePGMToImage<irtkVector3D<short> >;
template class irtkFilePGMToImage<irtkVector3D<float> >;
template class irtkFilePGMToImage<irtkVector3D<double> >;
