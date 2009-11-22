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

const char *irtkImageToFilePGM::NameOfClass()
{
  return "irtkImageToFilePGM";
}

void irtkImageToFilePGM::Initialize()
{
  unsigned int i;
  char header[255];

  // Initialize base class
  this->irtkImageToFile::Initialize();

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
  _start = i;
}

void irtkImageToFilePGM::Run()
{
  int i;
  double min, max;

  // Initialize filter
  this->Initialize();

  // Find dynamic range in image
  this->_input->GetMinMaxAsDouble(&min, &max);

  switch (this->_input->GetScalarType()) {
  case IRTK_VOXEL_CHAR: {
      // Copy voxels to buffer
      unsigned char *buffer = new unsigned char[this->_input->GetNumberOfVoxels()];

      char *ptr = (char *)this->_input->GetScalarPointer();
      for (i = 0; i < this->_input->GetNumberOfVoxels(); i++) {
        buffer[i] = round(255 * (*ptr - min) / (double)(max - min));
        ptr++;
      }

      // Write data
      this->WriteAsUChar(buffer, this->_input->GetNumberOfVoxels(), this->_start);

      // Destroy buffer
      delete []buffer;
      
      break;
    }
  case IRTK_VOXEL_UNSIGNED_CHAR: {
      // Write data
      this->WriteAsUChar((unsigned char *)this->_input->GetScalarPointer(), this->_input->GetNumberOfVoxels(), this->_start);
      
      break;
    }
  case IRTK_VOXEL_SHORT: {
      // Copy voxels to buffer
      unsigned char *buffer = new unsigned char[this->_input->GetNumberOfVoxels()];

      short *ptr = (short *)this->_input->GetScalarPointer();
      for (i = 0; i < this->_input->GetNumberOfVoxels(); i++) {
        buffer[i] = round(255 * (*ptr - min) / (double)(max - min));
        ptr++;
      }

      // Write data
      this->WriteAsUChar(buffer, this->_input->GetNumberOfVoxels(), this->_start);

      // Destroy buffer
      delete []buffer;
      
      break;
    }
  case IRTK_VOXEL_FLOAT: {
      // Copy voxels to buffer
      unsigned char *buffer = new unsigned char[this->_input->GetNumberOfVoxels()];

      float *ptr = (float *)this->_input->GetScalarPointer();
      for (i = 0; i < this->_input->GetNumberOfVoxels(); i++) {
        buffer[i] = round(255 * (*ptr - min) / (double)(max - min));
        ptr++;
      }

      // Write data
      this->WriteAsUChar(buffer, this->_input->GetNumberOfVoxels(), this->_start);

      // Destroy buffer
      delete []buffer;
      
      break;
    }
  case IRTK_VOXEL_DOUBLE: {
      // Copy voxels to buffer
      unsigned char *buffer = new unsigned char[this->_input->GetNumberOfVoxels()];

      double *ptr = (double *)this->_input->GetScalarPointer();
      for (i = 0; i < this->_input->GetNumberOfVoxels(); i++) {
        buffer[i] = round(255 * (*ptr - min) / (double)(max - min));
        ptr++;
      }

      // Write data
      this->WriteAsUChar(buffer, this->_input->GetNumberOfVoxels(), this->_start);

      // Destroy buffer
      delete []buffer;
      
      break;
    }
  default:
      cerr << "irtkImageToFilePGM::Run(): Unknown voxel type" << endl;
      exit(1);
  }

  // Finalize filter
  this->Finalize();
}



