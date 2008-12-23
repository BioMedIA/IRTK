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

#include <irtkImageToFile.h>

#include <irtkNIFTI.h>

template <class VoxelType> irtkImageToFileNIFTI<VoxelType>::irtkImageToFileNIFTI() : irtkImageToFile<VoxelType>()
{
  this->_headername = NULL;
}

template <class VoxelType> irtkImageToFileNIFTI<VoxelType>::~irtkImageToFileNIFTI()
{
  if (this->_headername != NULL) free(this->_headername);
}

template <class VoxelType> void irtkImageToFileNIFTI<VoxelType>::SetOutput(const char *name)
{
  int length;
  char *imagename;

  if (this->_headername != NULL) free(this->_headername);
  if (strstr(name, ".gz") == NULL) {
    this->_headername = strdup(name);
    imagename   = strdup(name);
    length      = strlen(name);
    imagename[length-1] = 'i';
    imagename[length-2] = 'i';
    imagename[length-3] = 'n';
  } else {
    this->_headername = strdup(name);
    imagename   = strdup(name);
    length      = strlen(name);
    imagename[length-4] = 'i';
    imagename[length-5] = 'i';
    imagename[length-6] = 'n';
  }
  this->irtkImageToFile<VoxelType>::SetOutput(imagename);
}

template <class VoxelType> const char *irtkImageToFileNIFTI<VoxelType>::NameOfClass()
{
  return "irtkImageToFileNIFTI";
}

template <class VoxelType> void irtkImageToFileNIFTI<VoxelType>::Finalize()
{
  // Empty method
  return;
}

template <> void irtkImageToFileNIFTI<irtkBytePixel>::Initialize()
{
  int x, y, z, t;
  double xsize, ysize, zsize, tsize;
  irtkMatrix i2w;

  // Get image info
  x = this->_input->GetX();
  y = this->_input->GetY();
  z = this->_input->GetZ();
  t = this->_input->GetT();

  this->_input->GetPixelSize(&xsize, &ysize, &zsize, &tsize);

  i2w = this->_input->GetImageToWorldMatrix();

  // Init header
  _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_UINT8, i2w);
}

template <> void irtkImageToFileNIFTI<irtkGreyPixel>::Initialize()
{
  int x, y, z, t;
  double xsize, ysize, zsize, tsize;
  irtkMatrix i2w;

  // Get image info
  x = this->_input->GetX();
  y = this->_input->GetY();
  z = this->_input->GetZ();
  t = this->_input->GetT();

  this->_input->GetPixelSize(&xsize, &ysize, &zsize, &tsize);

  i2w = this->_input->GetImageToWorldMatrix();

  // Init and write header
  _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_INT16, i2w);

}

template <> void irtkImageToFileNIFTI<irtkRealPixel>::Initialize()
{
  int x, y, z, t;
  double xsize, ysize, zsize, tsize;
  irtkMatrix i2w;

  // Get image info
  x = this->_input->GetX();
  y = this->_input->GetY();
  z = this->_input->GetZ();
  t = this->_input->GetT();

  this->_input->GetPixelSize(&xsize, &ysize, &zsize, &tsize);

  i2w = this->_input->GetImageToWorldMatrix();

  // Init header
  _hdr.Initialize(x, y, z, t, xsize, ysize, zsize, tsize, NIFTI_TYPE_FLOAT32, i2w);
}

template <> void irtkImageToFileNIFTI<irtkRGBPixel>::Initialize()
{
  cerr << "irtkImageToFileNIFTI<irtkRGBPixel>::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileNIFTI<irtkVector3D<char> >::Initialize()
{
  cerr << "irtkImageToFileNIFTI<irtkVector3D<char> >::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileNIFTI<irtkVector3D<short> >::Initialize()
{
  cerr << "irtkImageToFileNIFTI<irtkVector3D<short> >::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileNIFTI<irtkVector3D<float> >::Initialize()
{
  cerr << "irtkImageToFileNIFTI<irtkVector3D<float> >::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFileNIFTI<irtkVector3D<double> >::Initialize()
{
  cerr << "irtkImageToFileNIFTI<irtkVector3D<double> >::Run: Not supported" << endl;
  exit(1);
}

template <class VoxelType> void irtkImageToFileNIFTI<VoxelType>::Run()
{
  // Initialize filter
  this->Initialize();

  // Set filename in _hdr
  struct nifti_1_header nhdr = nifti_convert_nim2nhdr(_hdr.nim);
  _hdr.nim = nifti_convert_nhdr2nim(nhdr, this->_output); // This sets fname and iname
  _hdr.nim->iname_offset = 352;                           // Some nifti versions lose this on the way!

  /// Set data pointer in nifti image struct
  _hdr.nim->data = this->_input->GetPointerToVoxels();

  // Write hdr and data
  nifti_image_write(_hdr.nim);

  // Finalize filter
  this->Finalize();
}


template class irtkImageToFileNIFTI<irtkBytePixel>;
template class irtkImageToFileNIFTI<irtkGreyPixel>;
template class irtkImageToFileNIFTI<irtkRealPixel>;
template class irtkImageToFileNIFTI<irtkRGBPixel>;
template class irtkImageToFileNIFTI<irtkVector3D<char> >;
template class irtkImageToFileNIFTI<irtkVector3D<short> >;
template class irtkImageToFileNIFTI<irtkVector3D<float> >;
template class irtkImageToFileNIFTI<irtkVector3D<double> >;

#endif
