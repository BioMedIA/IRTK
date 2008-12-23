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

template <class VoxelType> irtkFileANALYZEToImage<VoxelType>::irtkFileANALYZEToImage()
{
  _headername = NULL;

  // Analyze-specific
  this->_reflectY = True;
}

template <class VoxelType> irtkFileANALYZEToImage<VoxelType>::~irtkFileANALYZEToImage()
{
  if (this->_headername != NULL) free(this->_headername);
}

template <class VoxelType> const char *irtkFileANALYZEToImage<VoxelType>::NameOfClass()
{
  return "irtkFileANALYZEToImage";
}

template <class VoxelType> int irtkFileANALYZEToImage<VoxelType>::CheckHeader(const char *filename)
{
  if ((strstr(basename2(filename), ".hdr") == NULL) && (strstr(basename2(filename), ".HDR") == NULL)) {
    return False;
  } else {
    return True;
  }
}

template <class VoxelType> void irtkFileANALYZEToImage<VoxelType>::SetInput(const char *headername)
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

  this->irtkFileToImage<VoxelType>::SetInput(imagename);
}

template <class VoxelType> void irtkFileANALYZEToImage<VoxelType>::ReadHeader()
{
  int i;
  irtkANALYZEHeader hdr;

  // Read header
  hdr.Read(this->_headername);

  // Copy header information
  if (hdr.dims[0] == 3) {
    this->_x     = hdr.dims[1];
    this->_y     = hdr.dims[2];
    this->_z     = hdr.dims[3];
    this->_t     = 1;
    this->_xsize = hdr.pixdims[1];
    this->_ysize = hdr.pixdims[2];
    this->_zsize = hdr.pixdims[3];
    this->_tsize = 1;
  } else {
    if (hdr.dims[0] == 4) {
      this->_x     = hdr.dims[1];
      this->_y     = hdr.dims[2];
      this->_z     = hdr.dims[3];
      this->_xsize = hdr.pixdims[1];
      this->_ysize = hdr.pixdims[2];
      this->_zsize = hdr.pixdims[3];
      if (hdr.dims[4] == 0) {
        this->_t     = 1;
        this->_tsize = 1;
      } else {
        this->_t     = hdr.dims[4];
        this->_tsize = hdr.pixdims[4];
      }
    } else {
      cerr << "irtkFileANALYZEToImage<VoxelType>::ReadHeader: Invalid no. of dimensions in image" << endl;
      exit(1);
    }
  }

  // If the size of header is not 348 we need to invert swapping
  if (hdr.sizeof_hdr != 348) {
    this->_swapped = !this->_swapped;
  }

  switch (hdr.data_type) {
  case ANALYZE_UNSIGNED_CHAR:
    this->_type  = VOXEL_U_CHAR;
    this->_bytes = 1;
    break;
  case ANALYZE_SIGNED_SHORT:
    this->_type  = VOXEL_SHORT;
    this->_bytes = 2;
    break;
  case ANALYZE_FLOAT:
    this->_type  = VOXEL_FLOAT;
    this->_bytes = 4;
    break;
  default:
    cerr << this->NameOfClass() << ": Data type " << hdr.data_type << " not supported, trying signed short data type" << endl;
    this->_type  = VOXEL_SHORT;
    this->_bytes = 2;
  }

  // Allocate memory for data address
  if (this->_addr != 0) delete []this->_addr;
  this->_addr = new int[this->_z];

  // Calculate data address
  for (i = 0; i < this->_z; i++) {
    this->_addr[i] = i * this->_x * this->_y * this->_bytes;
  }
}

template class irtkFileANALYZEToImage<irtkBytePixel>;
template class irtkFileANALYZEToImage<irtkGreyPixel>;
template class irtkFileANALYZEToImage<irtkRealPixel>;
template class irtkFileANALYZEToImage<irtkVector3D<char> >;
template class irtkFileANALYZEToImage<irtkVector3D<short> >;
template class irtkFileANALYZEToImage<irtkVector3D<float> >;
template class irtkFileANALYZEToImage<irtkVector3D<double> >;
