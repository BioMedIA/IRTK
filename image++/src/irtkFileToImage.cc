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

irtkFileToImage::irtkFileToImage()
{
  _type  = IRTK_VOXEL_UNKNOWN;
  _slope = 1.;
  _intercept = 0.;
  _reflectX = false;
  _reflectY = false;
  _reflectZ = false;
  _debug = true;
  _start = 0;
  _imagename = NULL;

}

irtkFileToImage::~irtkFileToImage()
{
  _type  = IRTK_VOXEL_UNKNOWN;
  _slope = 1.;
  _intercept = 0.;
  _reflectX = false;
  _reflectY = false;
  _reflectZ = false;
  _debug = true;
  _start = 0;
  if (_imagename != NULL) 
      _imagename = NULL;
}

irtkFileToImage *irtkFileToImage::New(const char *imagename)
{
  irtkFileToImage *reader = NULL;

  // Check format for GIPL
  if (irtkFileGIPLToImage::CheckHeader(imagename)) {
    reader = new irtkFileGIPLToImage;
    reader->SetInput(imagename);
    return reader;
  }

#ifdef HAS_NIFTI
  // Check format for NIFTI
  if (irtkFileNIFTIToImage::CheckHeader(imagename)) {
    reader = new irtkFileNIFTIToImage;
    reader->SetInput(imagename);
    return reader;
  }
#endif

  // Check format for ANALYZE
  if (irtkFileANALYZEToImage::CheckHeader(imagename)) {
    reader = new irtkFileANALYZEToImage;
    reader->SetInput(imagename);
    return reader;
  }

  // Check format for VTK
  if (irtkFileVTKToImage::CheckHeader(imagename)) {
    reader = new irtkFileVTKToImage;
    reader->SetInput(imagename);
    return reader;
  }

  // Check format for PGM
  if (irtkFilePGMToImage::CheckHeader(imagename)) {
    reader = new irtkFilePGMToImage;
    reader->SetInput(imagename);
    return reader;
  }

#ifdef HAS_OPENCV
  if (irtkFileOpenCVToImage::CheckHeader(imagename)) {
      reader = new irtkFileOpenCVToImage;
      reader->SetInput(imagename);
      return reader;
  }
#endif

  // Check for error
  if (reader == NULL) {
      stringstream msg;
      msg << "irtkFileToImage::New: Unknown file format " << imagename
          << endl;
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
  }

  return reader;
}

int irtkFileToImage::GetDebugFlag()
{
  return _debug;
}

void irtkFileToImage::PutDebugFlag(int debug)
{
  _debug = debug;
}

const char *irtkFileToImage::NameOfClass()
{
  return "irtkFileToImage";
}

void irtkFileToImage::SetInput(const char *imagename)
{
  // Close old file
  this->Close();

  // Copy new file name
  _imagename = imagename;

  // Open new file for reading
  this->Open(_imagename);

  // Read header
  this->ReadHeader();
}

irtkImage *irtkFileToImage::GetOutput()
{
  irtkImage *output = NULL;

  // Bring image to correct size
  switch (_type) {
  case IRTK_VOXEL_CHAR: {

      // Allocate image
      output = new irtkGenericImage<char>(_attr);

      // Read data
      this->ReadAsChar((char *)output->GetScalarPointer(), output->GetNumberOfVoxels(), _start);

    }
    break;

  case IRTK_VOXEL_UNSIGNED_CHAR: {

      // Allocate image
      output = new irtkGenericImage<unsigned char>(_attr);

      // Read data
      this->ReadAsUChar((unsigned char *)output->GetScalarPointer(), output->GetNumberOfVoxels(), _start);

    }
    break;

  case IRTK_VOXEL_SHORT: {

      // Allocate image
      output = new irtkGenericImage<short>(_attr);

      // Read data
      this->ReadAsShort((short *)output->GetScalarPointer(), output->GetNumberOfVoxels(), _start);

    }
    break;

  case IRTK_VOXEL_UNSIGNED_SHORT: {

      // Allocate image
      output = new irtkGenericImage<unsigned short>(_attr);

      // Read data
      this->ReadAsUShort((unsigned short *)output->GetScalarPointer(), output->GetNumberOfVoxels(), _start);

    }
    break;
    
  case IRTK_VOXEL_INT: {

      // Allocate image
      output = new irtkGenericImage<int>(_attr);

      // Read data
      this->ReadAsInt((int *)output->GetScalarPointer(), output->GetNumberOfVoxels(), _start);

    }
    break;

  case IRTK_VOXEL_UNSIGNED_INT: {
		  // Allocate image
		  output = new irtkGenericImage<unsigned int>(_attr);

		  // Read data
		  this->ReadAsUInt((unsigned int *)output->GetScalarPointer(), output->GetNumberOfVoxels(), _start);
  	  }
  	  break;
  	  
  case IRTK_VOXEL_FLOAT: {

      // Allocate image
      output = new irtkGenericImage<float>(_attr);

      // Read data
      this->ReadAsFloat((float *)output->GetScalarPointer(), output->GetNumberOfVoxels(), _start);

    }
    break;

  case IRTK_VOXEL_DOUBLE: {

      // Allocate image
      output = new irtkGenericImage<double>(_attr);

      // Read data
      this->ReadAsDouble((double *)output->GetScalarPointer(), output->GetNumberOfVoxels(), _start);

    }
    break;

   default:
      cout << "irtkFileToImage::GetOutput: Unknown voxel type" << endl;
  }

  // Reflect if necessary
  if (_reflectX == true) output->ReflectX();
  if (_reflectY == true) output->ReflectY();
  if (_reflectZ == true) output->ReflectZ();

  return output;
}

double irtkFileToImage::GetSlope()
{
	return this->_slope;
}

double irtkFileToImage::GetIntercept()
{
	return this->_intercept;
}

void irtkFileToImage::Print()
{
  cout << "Name of class is " << this->NameOfClass() << endl;
  cout << "File name is " << _imagename << endl;
  cout << "Image dimensions are " << _attr._x << " " << _attr._y << " " << _attr._z << " " << _attr._t << endl;
  cout << "Image has " << _bytes << " bytes per voxel" << endl;
  cout << "Voxel dimensions are " << _attr._dx << " " << _attr._dy << " "
  << _attr._dz << " " << _attr._dt << endl;
  cout << "Voxel type is ";
  switch (_type) {
    case IRTK_VOXEL_CHAR:
      cout << "char" << endl;
      break;
    case IRTK_VOXEL_UNSIGNED_CHAR:
      cout << "unsigned char" << endl;
      break;
    case IRTK_VOXEL_SHORT:
      cout << "short" << endl;
      break;
    case IRTK_VOXEL_UNSIGNED_SHORT:
      cout << "unsigned short" << endl;
      break;
    case IRTK_VOXEL_INT:
      cout << "int" << endl;
      break;
    case IRTK_VOXEL_UNSIGNED_INT:
      cout << "unsigned int" << endl;
      break;
    case IRTK_VOXEL_FLOAT:
      cout << "float" << endl;
      break;
    case IRTK_VOXEL_DOUBLE:
      cout << "double" << endl;
      break;
    default:
      cout << "unknown" << endl;
  }
}

void irtkFileToImage::Debug(char *message)
{
  if (_debug) cerr << message << endl;
}

int irtkFileToImage::GetDataType()
{
	return this->_type;
}

