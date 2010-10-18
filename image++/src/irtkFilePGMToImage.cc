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

const char *irtkFilePGMToImage::NameOfClass()
{
  return "irtkFilePGMToImage";
}

int irtkFilePGMToImage::CheckHeader(const char *filename)
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
    return false;
  } else {
    return true;
  }
}

void irtkFilePGMToImage::ReadHeader()
{
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
  sscanf(buffer, "%d %d", &this->_attr._x, &this->_attr._y);

  // Ignore maximum greyvalue, skip comments
  do {
    this->ReadAsString(buffer, 255);
  } while (buffer[0] == '#');

  // PGM files support only 2D images, so set z and t to 1
  this->_attr._z = 1;
  this->_attr._t = 1;

  // PGM files do not have voxel dimensions, so set them to default values
  this->_attr._dx = 1;
  this->_attr._dy = 1;
  this->_attr._dz = 1;
  this->_attr._dt = 1;

  // PGM files have voxels which are unsigned char
  this->_type  = IRTK_VOXEL_UNSIGNED_CHAR;
  this->_bytes = 1;

  // Data starts here
  this->_start = this->Tell();
}

