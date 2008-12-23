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

#include <irtkVTK.h>

template <class VoxelType> const char *irtkFileVTKToImage<VoxelType>::NameOfClass()
{
  return "irtkFileVTKToImage";
}

template <class VoxelType> int irtkFileVTKToImage<VoxelType>::CheckHeader(const char *filename)
{
  char buffer[255];

  // Create file stream
  irtkCifstream from;

  // Open new file for reading
  from.Open(filename);

  // Read header
  from.ReadAsString(buffer, 255);

  // Close file
  from.Close();

  // Check header
  if ((strcmp(buffer, VTK_MAGIC1) != 0) &&
      (strcmp(buffer, VTK_MAGIC2) != 0) &&
      (strcmp(buffer, VTK_MAGIC3) != 0)) {
    return False;
  } else {
    return True;
  }
}

template <class VoxelType> void irtkFileVTKToImage<VoxelType>::ReadHeader()
{
  int i;
  char type[255], dummy[255], buffer[255];

  // Initialize
  this->_type  = VOXEL_UNKNOWN;
  this->_bytes = 0;

  // Read header
  this->ReadAsString(buffer, 255);

  // Check header
  if ((strcmp(buffer, VTK_MAGIC1) != 0) &&
      (strcmp(buffer, VTK_MAGIC2) != 0) &&
      (strcmp(buffer, VTK_MAGIC3) != 0)) {
    cerr << this->NameOfClass() << "::Read_Header: Can't read magic number"
         << buffer << endl;
    exit(1);
  }

  // Skip comments
  this->ReadAsString(buffer, 255);

  // Read file type
  this->ReadAsString(buffer, 255);

  // Read dataset type
  this->ReadAsString(buffer, 255);

  // Check whether dataset is a structured points data set
  if (strcmp(buffer, "DATASET STRUCTURED_POINTS") != 0) {
    cerr << this->NameOfClass() << "::Read_Header: Data set is not in "
         << "structured points format - " << buffer << endl;
    exit(1);
  }

  // Read image dimensions type
  this->ReadAsString(buffer, 255);

  // Parse image dimensions
  sscanf(buffer, "%s %d %d %d", dummy, &this->_x, &this->_y, &this->_z);
  if (strcmp(dummy, "DIMENSIONS") != 0) {
    cerr << this->NameOfClass() << "::Read_Header: Can't find image dimensions"
         << endl;
  }

  // Read origin
  this->ReadAsString(buffer, 255);

  // Parse origin, but ignore
  sscanf(buffer, "%s %lf %lf %lf", dummy, &this->_xsize, &this->_ysize, &this->_zsize);
  if (strcmp(dummy, "ORIGIN") == 0) {

    // Read voxel dimensions
    this->ReadAsString(buffer, 255);

    // Parse voxel dimensions
    sscanf(buffer, "%s %lf %lf %lf", dummy, &this->_xsize, &this->_ysize, &this->_zsize);
    if ((strcmp(dummy, "SPACING") != 0) &&
        (strcmp(dummy, "ASPECT_RATIO") != 0)) {
      cerr << this->NameOfClass()
           << "::Read_Header: Can't find voxel dimensions" << endl;
    }
  } else {

    // Read voxel dimensions
    this->ReadAsString(buffer, 255);
  }

  // Read no. of points
  this->ReadAsString(buffer, 255);

  // Parse no. of points, but ignore
  sscanf(buffer, "%s %*d", dummy);
  if (strcmp(dummy, "POINT_DATA") != 0) {
    cerr << this->NameOfClass() << "::Read_Header: Can't find no. of points"
         << endl;
  }

  // Read voxel type
  this->ReadAsString(buffer, 255);

  // Parse voxel type
  sscanf(buffer, "%s %*s %s", dummy, type);

  if (strcmp(dummy, "SCALARS") != 0) {
    cerr << this->NameOfClass() << "::Read_Header: Can't find voxel type"
         << endl;
  }
  if (strcmp(type, VTK_DATA_CHAR) == 0) {
    this->_type  = VOXEL_CHAR;
    this->_bytes = 1;
  }
  if (strcmp(type, VTK_DATA_U_CHAR) == 0) {
    this->_type  = VOXEL_U_CHAR;
    this->_bytes = 1;
  }
  if (strcmp(type, VTK_DATA_SHORT) == 0) {
    this->_type  = VOXEL_SHORT;
    this->_bytes = 2;
  }
  if (strcmp(type, VTK_DATA_U_SHORT) == 0) {
    this->_type  = VOXEL_U_SHORT;
    this->_bytes = 2;
  }
  if (strcmp(type, VTK_DATA_FLOAT) == 0) {
    this->_type  = VOXEL_FLOAT;
    this->_bytes = 4;
  }
  if (this->_type == VOXEL_UNKNOWN) {
    cerr << this->NameOfClass() << "::Read_Header: Unknown voxel type"
         << endl;
  }

  // VTK only supports 3D images
  this->_t     = 1;
  this->_tsize = 1;

  // Read lookup table but ignore
  this->ReadAsString(buffer, 255);

  // Allocate memory for data address
  if (this->_addr != 0) delete []this->_addr;
  this->_addr = new int[this->_z];

  // Calculate data address
  for (i = 0; i < this->_z; i++) {
    this->_addr[i] = int(this->Tell()) + i * this->_x * this->_y * this->_bytes;
  }
}

template class irtkFileVTKToImage<irtkBytePixel>;
template class irtkFileVTKToImage<irtkGreyPixel>;
template class irtkFileVTKToImage<irtkRealPixel>;
template class irtkFileVTKToImage<irtkVector3D<char> >;
template class irtkFileVTKToImage<irtkVector3D<short> >;
template class irtkFileVTKToImage<irtkVector3D<float> >;
template class irtkFileVTKToImage<irtkVector3D<double> >;


