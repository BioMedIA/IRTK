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

template <class VoxelType> const char *irtkImageToFilePGM<VoxelType>::NameOfClass()
{
  return "irtkImageToFilePGM";
}

template <class VoxelType> void irtkImageToFilePGM<VoxelType>::Initialize()
{
  unsigned int i;
  char header[255];

  // Initialize base class
  this->irtkImageToFile<VoxelType>::Initialize();

  if (this->_input->GetZ() != 1) {
    cerr << this->NameOfClass() << " supports only 2D images" << endl;
    exit(1);
  }

  // Construct header
  sprintf(header, "%s\n# Created by %s\n%d %d\n255\n", PGM_MAGIC,
          this->NameOfClass(), this->_input->GetX(), this->_input->GetY());

  // Write header
  for (i = 0; i < strlen(header); i++) {
    this->WriteAsChar(header[i], i);
  }

  // Calculate data address
  this->_addr[0] = i;
}

template <> void irtkImageToFilePGM<irtkBytePixel>::Run()
{
  // Initialize filter
  this->Initialize();

  // Write data
  this->WriteAsUChar((irtkBytePixel *)this->_input->GetPointerToVoxels(),
                     this->_input->GetNumberOfVoxels(), this->_addr[0]);

  // Finalize filter
  this->Finalize();
}

template <> void irtkImageToFilePGM<irtkGreyPixel>::Run()
{
  int i;
  irtkGreyPixel *ptr, min, max;

  // Initialize filter
  this->Initialize();

  // Create buffer
  unsigned char *buffer = new unsigned char[this->_input->GetNumberOfVoxels()];

  // Find dynamic range in image
  this->_input->GetMinMax(&min, &max);

  // Copy voxels to buffer
  ptr = this->_input->GetPointerToVoxels();
  for (i = 0; i < this->_input->GetNumberOfVoxels(); i++) {
    buffer[i] = round(255 * (*ptr - min) / (double)(max - min));
    ptr++;
  }

  // Write data
  this->WriteAsUChar(buffer, this->_input->GetNumberOfVoxels(), this->_addr[0]);

  // Destroy buffer
  delete []buffer;

  // Finalize filter
  this->Finalize();
}

template <> void irtkImageToFilePGM<irtkRealPixel>::Run()
{
  int i;
  irtkRealPixel *ptr, min, max;

  // Initialize filter
  this->Initialize();

  // Create buffer
  unsigned char *buffer = new unsigned char[this->_input->GetNumberOfVoxels()];

  // Find dynamic range in image
  this->_input->GetMinMax(&min, &max);

  // Copy voxels to buffer
  ptr = this->_input->GetPointerToVoxels();
  for (i = 0; i < this->_input->GetNumberOfVoxels(); i++) {
    buffer[i] = round(255 * (*ptr - min) / (double)(max - min));
    ptr++;
  }

  // Write data
  this->WriteAsUChar(buffer, this->_input->GetNumberOfVoxels(), this->_addr[0]);

  // Destroy buffer
  delete []buffer;

  // Finalize filter
  this->Finalize();
}

template <> void irtkImageToFilePGM<irtkRGBPixel>::Run()
{
  cerr << "irtkImageToFilePGM<irtkRGBPixel>::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFilePGM<irtkVector3D<char> >::Run()
{
  cerr << "irtkImageToFilePGM<irtkVector3D<char> >::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFilePGM<irtkVector3D<short> >::Run()
{
  cerr << "irtkImageToFilePGM<irtkVector3D<short> >::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFilePGM<irtkVector3D<float> >::Run()
{
  cerr << "irtkImageToFilePGM<irtkVector3D<float> >::Run: Not supported" << endl;
  exit(1);
}

template <> void irtkImageToFilePGM<irtkVector3D<double> >::Run()
{
  cerr << "irtkImageToFilePGM<irtkVector3D<double> >::Run: Not supported" << endl;
  exit(1);
}

template class irtkImageToFilePGM<irtkBytePixel>;
template class irtkImageToFilePGM<irtkGreyPixel>;
template class irtkImageToFilePGM<irtkRealPixel>;
template class irtkImageToFilePGM<irtkRGBPixel>;
template class irtkImageToFilePGM<irtkVector3D<char> >;
template class irtkImageToFilePGM<irtkVector3D<short> >;
template class irtkImageToFilePGM<irtkVector3D<float> >;
template class irtkImageToFilePGM<irtkVector3D<double> >;

