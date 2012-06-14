/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifdef HAS_NIFTI

#include <irtkImage.h>

#include <irtkFileToImage.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <irtkNIFTI.h>

// little helpers (from fslio.h)
mat33 nifti_mat44_to_mat33(mat44 x)
{
  mat33 y;
  int i,j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      y.m[i][j] = x.m[i][j];
    }
  }
  return y;
}


irtkFileNIFTIToImage::irtkFileNIFTIToImage()
{
  _headername = NULL;
}

irtkFileNIFTIToImage::~irtkFileNIFTIToImage()
{
  if (this->_headername != NULL) free(this->_headername);
}

const char *irtkFileNIFTIToImage::NameOfClass()
{
  return "irtkFileNIFTIToImage";
}

int irtkFileNIFTIToImage::CheckHeader(const char *filename)
{
  char magic_number[4];

  // Create file stream
  irtkCifstream from;

  // Open new file for reading
  from.Open(filename);

  // Read magic no
  from.ReadAsChar(magic_number, 4, 344);

  // Close file
  from.Close();

  // Check magic no
  if ((strcmp(magic_number, "n+1") == 0) || (strcmp(magic_number, "ni1") == 0)) {
    return true;
  } else {
    return false;
  }
}

void irtkFileNIFTIToImage::SetInput(const char *filename)
{
  int length;
  char imagename[255], magic_number[4];
  struct stat buf;

  // Create file stream
  irtkCifstream from;

  // Open new file for reading
  from.Open(filename);

  // Read magic no
  from.ReadAsChar(magic_number, 4, 344);

  // Close file
  from.Close();

  // Delete old file name
  if (this->_headername != NULL) free(this->_headername);

  // Copy new file name
  this->_headername = strdup(filename);

  // Check magic no
  if (strcmp(magic_number, "n+1") == 0) {
    this->irtkFileToImage::SetInput(filename);
    return;
  } else {
    if (strcmp(magic_number, "ni1") == 0) {
      // Parse img file name
      if (strstr(basename2(filename), ".hdr.gz") != NULL) {
        length = strlen(filename);
        sprintf(imagename, "%s", filename);
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
        length = strlen(filename);
        sprintf(imagename, "%s", filename);
        imagename[length-1] = 'g';
        imagename[length-2] = 'm';
        imagename[length-3] = 'i';

        // Check if uncompressed file exists
        if (stat(imagename, &buf) != 0) {
          sprintf(imagename, "%s.gz", filename);
          imagename[length-1] = 'g';
          imagename[length-2] = 'm';
          imagename[length-3] = 'i';

          // Check if gzip compressed file exists
          if (stat(imagename, &buf) != 0) {
            sprintf(imagename, "%s.Z", filename);
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
    } else {
      cerr << "irtkFileNIFTIToImage::SetInput: File format is not NIFTI" << endl;
      exit(1);
    }
  }
  this->irtkFileToImage::SetInput(imagename);
}

void irtkFileNIFTIToImage::ReadHeader()
{
  int i, order;
  float det;
  mat44 mat_44, mat_inv_44;
  mat33 mat_33;
  irtkNIFTIHeader hdr;

  // Read header
  hdr.Read(this->_headername);
#ifdef DEBUG
  hdr.Print();
#endif

  // Check dimension
  if (hdr.nim->dim[0] > 4) {
    cerr << "irtkFileNIFTIToImage::ReadHeader: Number of dimensions > 4 (Number of dimensions = " << hdr.nim->dim[0] << ") \n";
    exit(1);
  }

  // Check intent code
  /*
  if (hdr.nim->intent_code != 0) {
    cerr << "irtkFileNIFTIToImage::ReadHeader: Unknown intent_code = " <<
         hdr.nim->intent_code << endl;
    exit(1);
  }
  */

  // Check for swapping
#ifndef WORDS_BIGENDIAN
  if (hdr.nim->byteorder == 1) { //LSB_FIRST) {
    this->_swapped = !this->_swapped;
  }
#else
  if (hdr.nim->byteorder == MSB_FIRST) {
    this->_swapped = !this->_swapped;
  }
#endif

  // Check data scaling
  if (hdr.nim->scl_slope != 0) {
    this->_slope = hdr.nim->scl_slope;
  } else {
    this->_slope = 1.0;
  }
  this->_intercept = hdr.nim->scl_inter;

  // Copy header information
  // (from now on force all vox dims to be positive - LR info is in sform)
  this->_attr._x  = hdr.nim->dim[1];
  this->_attr._y  = hdr.nim->dim[2];
  this->_attr._z  = hdr.nim->dim[3];
  this->_attr._dx = fabs(hdr.nim->pixdim[1]);
  this->_attr._dy = fabs(hdr.nim->pixdim[2]);
  this->_attr._dz = fabs(hdr.nim->pixdim[3]);
  if (hdr.nim->dim[0] == 4) {
    this->_attr._t  = hdr.nim->dim[4];
    this->_attr._dt = fabs(hdr.nim->pixdim[4]);
  } else {
    this->_attr._t  = 1;
    this->_attr._dt = 1;
  }

  // Check which coordinate system to use
  if (hdr.nim->qform_code > 0) {
    // Access qform
    mat_44 = hdr.nim->qto_xyz;
    mat_inv_44 = hdr.nim->qto_ijk;
  } else if (hdr.nim->sform_code > 0) {
    mat_44     = hdr.nim->sto_xyz;
    mat_inv_44 = hdr.nim->sto_ijk;
  } else {

    // What's below should be in hdr.nim->qto_xyz

    // grid spacings along diagonal (default radiological convention)
    mat_44.m[0][0] = -fabs(hdr.nim->pixdim[1]);
    mat_44.m[1][1] =  fabs(hdr.nim->pixdim[2]);
    mat_44.m[2][2] =  fabs(hdr.nim->pixdim[3]);

    // off diagonal is zero
    mat_44.m[0][1] = mat_44.m[0][2] = 0.0;
    mat_44.m[1][0] = mat_44.m[1][2] = 0.0;
    mat_44.m[2][0] = mat_44.m[2][1] = 0.0;

    // Apart from (inverted) offset in last column for origin to become (0,0,0) below
    mat_44.m[0][3] =  this->_attr._dx*(hdr.nim->dim[1] - 1) / 2.0;
    mat_44.m[1][3] = -this->_attr._dy*(hdr.nim->dim[2] - 1) / 2.0;
    mat_44.m[2][3] = -this->_attr._dz*(hdr.nim->dim[3] - 1) / 2.0;

    // last row is always [ 0 0 0 1 ]
    mat_44.m[3][0] = mat_44.m[3][1] = mat_44.m[3][2] = 0.0; mat_44.m[3][3 ]= 1.0 ;

    // Invert
    mat_inv_44 = nifti_mat44_inverse(mat_44);
    //mat_inv_44 = mat44_inverse(mat_44);
  }

  // Get left/right order (no check for inconsistency between qform/sform since only one is used)
  mat_33 = nifti_mat44_to_mat33(mat_44);
  det = nifti_mat33_determ(mat_33);
  // det = mat33_determ(mat_33);
  if (det < 0.0) {
    order = NIFTI_RADIOLOGICAL;
  } else {
    order = NIFTI_NEUROLOGICAL;
  }

  // Set axis orientation, including zaxis.
  // Need to preserve sign of axis, hence use absolute pixel sizes for descaling.
  for (i = 0; i < 3; i++) {
    this->_attr._xaxis[i] = mat_44.m[i][0] / this->_attr._dx;
    this->_attr._yaxis[i] = mat_44.m[i][1] / this->_attr._dy;
    this->_attr._zaxis[i] = mat_44.m[i][2] / this->_attr._dz;
  }

  // Convert between nifti and irtk coordinate systems
  // See https://www.fmrib.ox.ac.uk/ibim/uploads/coordtransforms.pdf
  irtkMatrix D(4, 4), D_inv(4, 4), M(4, 4), R;
  for (int j = 0; j < 4; j++) {
    M(j, j) = 1;
    for (i = 0; i < 4; i++) {
      D(i, j)     = mat_44.m[i][j];
      D_inv(i, j) = mat_inv_44.m[i][j];
    }
  }
  for (i = 0; i < 3; i++) M(i, 3) = (hdr.nim->dim[i+1] - 1) / 2.0;
  R = D * M * D_inv;

  // Set image origin by adding q/sform offset to third column of R:
  this->_attr._xorigin  = R(0, 3) + mat_44.m[0][3];
  this->_attr._yorigin  = R(1, 3) + mat_44.m[1][3];
  this->_attr._zorigin  = R(2, 3) + mat_44.m[2][3];

  // Set data type and number of bytes per voxels
  switch (hdr.nim->datatype) {
  case NIFTI_TYPE_UINT8:
    this->_type  = IRTK_VOXEL_UNSIGNED_CHAR;
    this->_bytes = 1;
    break;
  case NIFTI_TYPE_INT16:
    this->_type  = IRTK_VOXEL_SHORT;
    this->_bytes = 2;
    break;
  case NIFTI_TYPE_UINT16:
    this->_type  = IRTK_VOXEL_UNSIGNED_SHORT;
    this->_bytes = 2;
    break;
  case NIFTI_TYPE_FLOAT32:
    this->_type  = IRTK_VOXEL_FLOAT;
    this->_bytes = 4;
    break;
  case NIFTI_TYPE_INT32:
    this->_type  = IRTK_VOXEL_INT;
    this->_bytes = 4;
    break;
  case NIFTI_TYPE_UINT32:
    this->_type  = IRTK_VOXEL_UNSIGNED_INT;
    this->_bytes = 4;
    break;
  case NIFTI_TYPE_FLOAT64:
    this->_type = IRTK_VOXEL_DOUBLE;
    this->_bytes = 8;
    break;
  default:
    cerr << this->NameOfClass() << ": Data type " << hdr.nim->datatype << " not supported, trying signed short data type" << endl;
    this->_type  = IRTK_VOXEL_SHORT;
    this->_bytes = 2;
  }

  // Data starts here 
  this->_start = hdr.nim->iname_offset;
}

void irtkFileNIFTIToImage::Print()
{

  int i;
  irtkMatrix mat(3, 3);
  irtkNIFTIHeader hdr;

  hdr.Read(this->_headername);

  for (i = 0; i < 3; i++) {
    mat(0, i) = this->_attr._xaxis[i] * this->_attr._dx;
    mat(1, i) = this->_attr._yaxis[i] * this->_attr._dy;
    mat(2, i) = this->_attr._zaxis[i] * this->_attr._dz;
  }

  cout << "Name of class is " << this->NameOfClass() << endl;

  // Check which coordinate system is used.
  cout << "Transformation specified in NIFTI header using ";
  if (hdr.nim->qform_code > 0) {
    cout << "quaternion." << endl;
  } else if (hdr.nim->sform_code > 0) {
    cout << "affine transformation." << endl;
  } else {
    cout << "default transformation." << endl;
  }

  if (mat.Det() < 0.0) {
    cout << "Order is radiological\n";
  } else {
    cout << "Order is neurological\n";
  }

  cout << "Scale slope = " << this->_slope << endl;
  cout << "Intercept   = " << this->_intercept << endl;

  this->irtkFileToImage::Print();

}

#endif
