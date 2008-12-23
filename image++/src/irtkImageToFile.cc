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

template <class VoxelType> irtkImageToFile<VoxelType>::irtkImageToFile()
{
  _input  = NULL;
  _output = NULL;
  _addr   = NULL;
  _reflectX = False;
  _reflectY = False;
  _reflectZ = False;
}

template <class VoxelType> irtkImageToFile<VoxelType>::~irtkImageToFile()
{
  _input  = NULL;
  if (_output != NULL) free(_output);
  _output = NULL;
  if (_addr != NULL) delete []_addr;
  _addr   = NULL;
}

template <> irtkImageToFile<irtkRGBPixel> *irtkImageToFile<irtkRGBPixel>::New(const char *imagename)
{
  irtkImageToFile *writer = NULL;

#ifdef HAS_PNG
  writer = new irtkImageToFilePNG<irtkRGBPixel>;
  writer->SetOutput(imagename);
#else
  cerr << "irtkImageToFile<irtkRGBPixel>::New: No support for PNG file format compiled." << endl;
  exit(1);
#endif

  return writer;
}

template <class VoxelType> irtkImageToFile<VoxelType> *irtkImageToFile<VoxelType>::New(const char *imagename)
{
  irtkImageToFile *writer = NULL;

  // Check format for GIPL
  if (strstr(imagename, ".gipl") != NULL) {
    writer = new irtkImageToFileGIPL<VoxelType>;
    writer->SetOutput(imagename);
  }

  // Check format for VTK
  if (strstr(imagename, ".vtk") != NULL) {
    writer = new irtkImageToFileVTK<VoxelType>;
    writer->SetOutput(imagename);
  }

  // Check format for PGM
  if (strstr(imagename, ".pgm") != NULL) {
    writer = new irtkImageToFilePGM<VoxelType>;
    writer->SetOutput(imagename);
  }

  // Check format for ANALYZE
  if (strstr(imagename, ".hdr") != NULL) {
    writer = new irtkImageToFileANALYZE<VoxelType>;
    writer->SetOutput(imagename);
  }

#ifdef HAS_NIFTI
  // Check format for NIFTI
  if (strstr(imagename, ".nii") != NULL) {
    writer = new irtkImageToFileNIFTI<VoxelType>;
    writer->SetOutput(imagename);
  }
#endif

  // Check for default file format
  if (writer == NULL) {
    writer = new irtkImageToFileGIPL<VoxelType>;
    writer->SetOutput(imagename);
  }

  return writer;
}

template <class VoxelType> const char *irtkImageToFile<VoxelType>::NameOfClass()
{
  return "irtkImageToFile";
}

template <class VoxelType> void irtkImageToFile<VoxelType>::SetInput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << this->NameOfClass() << "::SetInput: Output is not an image\n";
    exit(1);
  }
}

template <class VoxelType> void irtkImageToFile<VoxelType>::SetOutput(const char *name)
{
  if (name != NULL) {
    _output = strdup(name);
  } else {
    cerr << this->NameOfClass() << "::SetInput: Output is not a filename\n";
    exit(1);
  }
}

template <class VoxelType> void irtkImageToFile<VoxelType>::Initialize()
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

  if (_addr != NULL) delete []_addr;
  _addr = new int[_input->GetZ()];

  // Reflect if necessary
  if (_reflectX == True) _input->ReflectX();
  if (_reflectY == True) _input->ReflectY();
  if (_reflectZ == True) _input->ReflectZ();
}

template <class VoxelType> void irtkImageToFile<VoxelType>::Finalize()
{
  // Close file
  this->Close();

  // Reflect back if necessary
  if (_reflectX == True) _input->ReflectX();
  if (_reflectY == True) _input->ReflectY();
  if (_reflectZ == True) _input->ReflectZ();
}

template <> void irtkImageToFile<irtkBytePixel>::Run()
{
  // Initialize filter
  this->Initialize();

  // Write data
  this->WriteAsUChar((irtkBytePixel *)_input->GetPointerToVoxels(),
                     _input->GetNumberOfVoxels(), _addr[0]);

  // Finalize filter
  this->Finalize();
}

template <> void irtkImageToFile<irtkGreyPixel>::Run()
{
  // Initialize filter
  this->Initialize();

  // Write data
  this->WriteAsShort((irtkGreyPixel *)_input->GetPointerToVoxels(),
                     _input->GetNumberOfVoxels(), _addr[0]);

  // Finalize filter
  this->Finalize();
}

template <> void irtkImageToFile<irtkRealPixel>::Run()
{
  // Initialize filter
  this->Initialize();

  // Write data
  this->WriteAsFloat((irtkRealPixel *)_input->GetPointerToVoxels(),
                     _input->GetNumberOfVoxels(), _addr[0]);

  // Finalize filter
  this->Finalize();
}

template <> void irtkImageToFile<irtkRGBPixel>::Run()
{
}

template <> void irtkImageToFile<irtkVector3D<char> >::Run()
{
  // Initialize filter.
  this->Initialize();

  // Write data.
  unsigned long index = 0;
  irtkVector3D<char>* pVoxels = _input->GetPointerToVoxels();

  for (int z = 0; z < _input->GetZ(); z++) {
    for (int y = 0; y < _input->GetY(); y++) {
      for (int x = 0; x < _input->GetX(); x++) {
        this->WriteAsChar(&pVoxels[index]._x, 1);
        this->WriteAsChar(&pVoxels[index]._y, 1);
        this->WriteAsChar(&pVoxels[index]._z, 1);

        index++;
      }
    }
  }

  // Finalize filter.
  this->Finalize();
}

template <> void irtkImageToFile<irtkVector3D<short> >::Run()
{
  // Initialize filter.
  this->Initialize();

  // Write data.
  unsigned long index = 0;
  irtkVector3D<short>* pVoxels = _input->GetPointerToVoxels();

  for (int z = 0; z < _input->GetZ(); z++) {
    for (int y = 0; y < _input->GetY(); y++) {
      for (int x = 0; x < _input->GetX(); x++) {
        this->WriteAsShort(&pVoxels[index]._x, 1);
        this->WriteAsShort(&pVoxels[index]._y, 1);
        this->WriteAsShort(&pVoxels[index]._z, 1);

        index++;
      }
    }
  }

  // Finalize filter.
  this->Finalize();
}

template <> void irtkImageToFile<irtkVector3D<float> >::Run()
{
  // Initialize filter.
  this->Initialize();

  // Write data.
  unsigned long index = 0;
  irtkVector3D<float>* pVoxels = _input->GetPointerToVoxels();

  for (int z = 0; z < _input->GetZ(); z++) {
    for (int y = 0; y < _input->GetY(); y++) {
      for (int x = 0; x < _input->GetX(); x++) {
        this->WriteAsFloat(&pVoxels[index]._x, 1);
        this->WriteAsFloat(&pVoxels[index]._y, 1);
        this->WriteAsFloat(&pVoxels[index]._z, 1);

        index++;
      }
    }
  }

  // Finalize filter.
  this->Finalize();
}

template <> void irtkImageToFile<irtkVector3D<double> >::Run()
{
  // Initialize filter.
  this->Initialize();

  // Write data.
  unsigned long index = 0;
  irtkVector3D<double>* pVoxels = _input->GetPointerToVoxels();

  for (int z = 0; z < _input->GetZ(); z++) {
    for (int y = 0; y < _input->GetY(); y++) {
      for (int x = 0; x < _input->GetX(); x++) {
        this->WriteAsDouble(&pVoxels[index]._x, 1);
        this->WriteAsDouble(&pVoxels[index]._y, 1);
        this->WriteAsDouble(&pVoxels[index]._z, 1);

        index++;
      }
    }
  }

  // Finalize filter.
  this->Finalize();
}

template class irtkImageToFile<irtkBytePixel>;
template class irtkImageToFile<irtkGreyPixel>;
template class irtkImageToFile<irtkRealPixel>;
template class irtkImageToFile<irtkRGBPixel>;
template class irtkImageToFile<irtkVector3D<char> >;
template class irtkImageToFile<irtkVector3D<short> >;
template class irtkImageToFile<irtkVector3D<float> >;
template class irtkImageToFile<irtkVector3D<double> >;

