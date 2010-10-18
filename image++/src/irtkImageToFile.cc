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

irtkImageToFile::irtkImageToFile()
{
  _input  = NULL;
  _output = NULL;
  _start   = 0;
  _reflectX = false;
  _reflectY = false;
  _reflectZ = false;
}

irtkImageToFile::~irtkImageToFile()
{
  _input  = NULL;
  if (_output != NULL) free(_output);
  _output = NULL;
  _start = 0;
}

irtkImageToFile *irtkImageToFile::New(const char *imagename)
{
  irtkImageToFile *writer = NULL;

#ifdef HAS_PNG

  // Check format for PNG
  if (strstr(imagename, ".png") != NULL) {
    writer = new irtkImageToFilePNG;
    writer->SetOutput(imagename);
  }

#endif
  
  // Check format for GIPL
  if (strstr(imagename, ".gipl") != NULL) {
    writer = new irtkImageToFileGIPL;
    writer->SetOutput(imagename);
  }

  // Check format for VTK
  if (strstr(imagename, ".vtk") != NULL) {
    writer = new irtkImageToFileVTK;
    writer->SetOutput(imagename);
  }

  // Check format for PGM
  if (strstr(imagename, ".pgm") != NULL) {
    writer = new irtkImageToFilePGM;
    writer->SetOutput(imagename);
  }

  // Check format for ANALYZE
  if (strstr(imagename, ".hdr") != NULL) {
    writer = new irtkImageToFileANALYZE;
    writer->SetOutput(imagename);
  }

#ifdef HAS_NIFTI
  // Check format for NIFTI
  if (strstr(imagename, ".nii") != NULL) {
    writer = new irtkImageToFileNIFTI;
    writer->SetOutput(imagename);
  }
#endif

  // Check for default file format
  if (writer == NULL) {
    writer = new irtkImageToFileGIPL;
    writer->SetOutput(imagename);
  }

  return writer;
}

const char *irtkImageToFile::NameOfClass()
{
  return "irtkImageToFile";
}

void irtkImageToFile::SetInput(irtkImage *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << this->NameOfClass() << "::SetInput: Output is not an image\n";
    exit(1);
  }
}

void irtkImageToFile::SetOutput(const char *name)
{
  if (name != NULL) {
    _output = strdup(name);
  } else {
    cerr << this->NameOfClass() << "::SetInput: Output is not a filename\n";
    exit(1);
  }
}

void irtkImageToFile::Initialize()
{
  // Check inputs and outputs
  if (_input == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input" << endl;
    exit(1);
  }
  if (_output == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no output" << endl;
    exit(1);
  }

  // Open file for writing
  this->Open(_output);

  // Reflect if necessary
  if (_reflectX == true) _input->ReflectX();
  if (_reflectY == true) _input->ReflectY();
  if (_reflectZ == true) _input->ReflectZ();
}

void irtkImageToFile::Finalize()
{
  // Close file
  this->Close();

  // Reflect back if necessary
  if (_reflectX == true) _input->ReflectX();
  if (_reflectY == true) _input->ReflectY();
  if (_reflectZ == true) _input->ReflectZ();
}

void irtkImageToFile::Run()
{
  // Initialize filter
  this->Initialize();

  // Write data
  switch (this->_input->GetScalarType()) {
  case IRTK_VOXEL_CHAR: {
      this->WriteAsChar((char *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_UNSIGNED_CHAR: {
      this->WriteAsUChar((unsigned char *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_SHORT: {
      this->WriteAsShort((short *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_UNSIGNED_SHORT: {
      this->WriteAsUShort((unsigned short *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_FLOAT: {
      this->WriteAsFloat((float *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_DOUBLE: {
      this->WriteAsDouble((double *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  default:
  	cerr << "irtkImageToFile::Run(): Unknown voxel type" << endl;
  	exit(1);
  }
  
  // Finalize filter
  this->Finalize();
}
