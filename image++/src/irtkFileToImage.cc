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

template <class VoxelType> irtkFileToImage<VoxelType>::irtkFileToImage()
{
  _x = 0;
  _y = 0;
  _z = 0;
  _t = 0;
  _xsize = 0;
  _ysize = 0;
  _zsize = 0;
  _tsize = 0;
  _xorigin = 0;
  _yorigin = 0;
  _zorigin = 0;
  _torigin = 0;
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;
  _type  = VOXEL_UNKNOWN;
  _slope = 1.;
  _intercept = 0.;
  _reflectX = False;
  _reflectY = False;
  _reflectZ = False;
  _debug = True;
  _addr  = NULL;
  _imagename = NULL;

}

template <class VoxelType> irtkFileToImage<VoxelType>::~irtkFileToImage()
{
  _x = 0;
  _y = 0;
  _z = 0;
  _t = 0;
  _xsize = 0;
  _ysize = 0;
  _zsize = 0;
  _tsize = 0;
  _xorigin = 0;
  _yorigin = 0;
  _zorigin = 0;
  _torigin = 0;
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;
  _bytes = 0;
  _type  = VOXEL_UNKNOWN;
  _debug = True;
  if (_addr != NULL) delete []_addr;
  if (_imagename != NULL) free(_imagename);
}

template <class VoxelType> irtkFileToImage<VoxelType> *irtkFileToImage<VoxelType>::New(const char *imagename)
{
  irtkFileToImage *reader = NULL;

  // Check format for GIPL
  if (irtkFileGIPLToImage<VoxelType>::CheckHeader(imagename)) {
    reader = new irtkFileGIPLToImage<VoxelType>;
    reader->SetInput(imagename);
    return reader;
  }

#ifdef HAS_NIFTI
  // Check format for NIFTI
  if (irtkFileNIFTIToImage<VoxelType>::CheckHeader(imagename)) {
    reader = new irtkFileNIFTIToImage<VoxelType>;
    reader->SetInput(imagename);
    return reader;
  }
#endif

  // Check format for ANALYZE
  if (irtkFileANALYZEToImage<VoxelType>::CheckHeader(imagename)) {
    reader = new irtkFileANALYZEToImage<VoxelType>;
    reader->SetInput(imagename);
    return reader;
  }

  // Check format for VTK
  if (irtkFileVTKToImage<VoxelType>::CheckHeader(imagename)) {
    reader = new irtkFileVTKToImage<VoxelType>;
    reader->SetInput(imagename);
    return reader;
  }

  // Check format for PGM
  if (irtkFilePGMToImage<VoxelType>::CheckHeader(imagename)) {
    reader = new irtkFilePGMToImage<VoxelType>;
    reader->SetInput(imagename);
    return reader;
  }

  // Check for error
  if (reader == NULL) {
    cerr << "irtkFileToImage::New: Unknown file format " << imagename
         << endl;
    exit(1);
  }

  return reader;
}

template <class VoxelType> int irtkFileToImage<VoxelType>::GetDebugFlag()
{
  return _debug;
}

template <class VoxelType> void irtkFileToImage<VoxelType>::PutDebugFlag(int debug)
{
  _debug = debug;
}

template <class VoxelType> const char *irtkFileToImage<VoxelType>::NameOfClass()
{
  return "irtkFileToImage";
}

template <class VoxelType> void irtkFileToImage<VoxelType>::SetInput(const char *imagename)
{
  // Close old file
  this->Close();

  // Delete old file name
  if (_imagename != NULL) free(_imagename);

  // Copy new file name
  _imagename = strdup(imagename);

  // Open new file for reading
  this->Open(_imagename);

  // Read header
  this->ReadHeader();
}

template <class VoxelType> void irtkFileToImage<VoxelType>::SetOutput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _output = image;
  } else {
    cerr << this->NameOfClass() << "::SetOutput: Output is not an image\n";
    exit(1);
  }
}

template <class VoxelType> void irtkFileToImage<VoxelType>::Run()
{
  char *data;
  int x, y, z, t;
  VoxelType *ptr1;

  // Bring image to correct size
  _output->Initialize(_x, _y, _z, _t, _xsize, _ysize, _zsize, _tsize, _xaxis, _yaxis, _zaxis);
  _output->PutOrigin(_xorigin, _yorigin, _zorigin, _torigin);

  // Allocate memory for tempory data
  data = new char[_x*_y*_z*_t*_bytes];

  // Get pointer to image data
  ptr1 = _output->GetPointerToVoxels();

  // Convert data
  switch (_type) {
  case VOXEL_CHAR: {
    char *ptr2 = &(data[0]);

    // Read data
    this->ReadAsChar((char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            *ptr1 = (VoxelType)(*ptr2 * _slope + _intercept);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_CHAR: {
    unsigned char *ptr2 = (unsigned char *)&(data[0]);

    // Read data
    this->ReadAsUChar((unsigned char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            *ptr1 = (VoxelType)(*ptr2 * _slope + _intercept);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data
    this->ReadAsShort((short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            *ptr1 = (VoxelType)(*ptr2 * _slope + _intercept);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_SHORT: {
    unsigned short *ptr2 = (unsigned short *)&(data[0]);

    // Read data
    this->ReadAsUShort((unsigned short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            *ptr1 = (VoxelType)(*ptr2 * _slope + _intercept);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_INT: {
    int *ptr2 = (int *)&(data[0]);

    // Read data
    this->ReadAsInt((int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            *ptr1 = (VoxelType)(*ptr2 * _slope + _intercept);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_INT: {
    unsigned int *ptr2 = (unsigned int *)&(data[0]);

    // Read data
    this->ReadAsUInt((unsigned int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            *ptr1 = (VoxelType)(*ptr2 * _slope + _intercept);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data
    this->ReadAsFloat((float *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            *ptr1 = (VoxelType)(*ptr2 * _slope + _intercept);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data
    this->ReadAsDouble((double *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            *ptr1 = (VoxelType)(*ptr2 * _slope + _intercept);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_CHAR: {
    char *ptr2 = (char *)&(data[0]);

    // Read data.
    this->ReadAsChar((char *)data, _x*_y*_z*3*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            double dx = *ptr2;
            double dy = *(ptr2 + 1);
            double dz = *(ptr2 + 2);

            *ptr1 = (VoxelType)(sqrt(dx*dx + dy*dy + dz*dz));
            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data.
    this->ReadAsShort((short *)data, _x*_y*_z*3*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            double dx = *ptr2;
            double dy = *(ptr2 + 1);
            double dz = *(ptr2 + 2);

            *ptr1 = (VoxelType)(sqrt(dx*dx + dy*dy + dz*dz));
            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data.
    this->ReadAsFloat((float *)data, _x*_y*_z*3*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            double dx = *ptr2;
            double dy = *(ptr2 + 1);
            double dz = *(ptr2 + 2);

            *ptr1 = (VoxelType)(sqrt(dx*dx + dy*dy + dz*dz));
            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data.
    this->ReadAsDouble((double *)data, _x*_y*_z*3*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            double dx = *ptr2;
            double dy = *(ptr2 + 1);
            double dz = *(ptr2 + 2);

            *ptr1 = (VoxelType)(sqrt(dx*dx + dy*dy + dz*dz));
            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  default:
    cout << "unknown" << endl;
  }

  // Delete data
  delete []data;

  // Reflect if necessary
  if (_reflectX == True) _output->ReflectX();
  if (_reflectY == True) _output->ReflectY();
  if (_reflectZ == True) _output->ReflectZ();
}

template <> void irtkFileToImage<irtkVector3D<char> >::Run()
{
  char *data;
  int x, y, z, t;
  irtkVector3D<char> *ptr1;

  // Bring image to correct size
  _output->Initialize(_x, _y, _z, _t, _xsize, _ysize, _zsize, _tsize, _xaxis, _yaxis, _zaxis);
  _output->PutOrigin(_xorigin, _yorigin, _zorigin, _torigin);

  // Allocate memory for tempory data
  data = new char[_x*_y*_z*_t*_bytes];

  // Get pointer to image data
  ptr1 = _output->GetPointerToVoxels();

  // Convert data
  switch (_type) {
  case VOXEL_CHAR: {
    char *ptr2 = &(data[0]);

    // Read data
    this->ReadAsChar((char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<char>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_CHAR: {
    unsigned char *ptr2 = (unsigned char *)&(data[0]);

    // Read data
    this->ReadAsUChar((unsigned char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<char>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data
    this->ReadAsShort((short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<char>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_SHORT: {
    unsigned short *ptr2 = (unsigned short *)&(data[0]);

    // Read data
    this->ReadAsUShort((unsigned short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<char>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_INT: {
    int *ptr2 = (int *)&(data[0]);

    // Read data
    this->ReadAsInt((int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<char>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_INT: {
    unsigned int *ptr2 = (unsigned int *)&(data[0]);

    // Read data
    this->ReadAsUInt((unsigned int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = *ptr2;
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data
    this->ReadAsFloat((float *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<char>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data
    this->ReadAsDouble((double *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<char>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_CHAR: {
    char *ptr2 = (char *)&(data[0]);

    // Read data.
    this->ReadAsChar((char *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<char>(*ptr2);
            ptr1->_y = static_cast<char>(*(ptr2 + 1));
            ptr1->_z = static_cast<char>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data.
    this->ReadAsShort((short *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<char>(*ptr2);
            ptr1->_y = static_cast<char>(*(ptr2 + 1));
            ptr1->_z = static_cast<char>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data.
    this->ReadAsFloat((float *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<char>(*ptr2);
            ptr1->_y = static_cast<char>(*(ptr2 + 1));
            ptr1->_z = static_cast<char>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data.
    this->ReadAsDouble((double *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<char>(*ptr2);
            ptr1->_y = static_cast<char>(*(ptr2 + 1));
            ptr1->_z = static_cast<char>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  default:
    cout << "unknown" << endl;
  }

  // Delete data
  delete []data;

  // Reflect if necessary
  if (_reflectX == True) _output->ReflectX();
  if (_reflectY == True) _output->ReflectY();
  if (_reflectZ == True) _output->ReflectZ();
}

template <> void irtkFileToImage<irtkVector3D<short> >::Run()
{
  char *data;
  int x, y, z, t;
  irtkVector3D<short> *ptr1;

  // Bring image to correct size
  _output->Initialize(_x, _y, _z, _t, _xsize, _ysize, _zsize, _tsize, _xaxis, _yaxis, _zaxis);
  _output->PutOrigin(_xorigin, _yorigin, _zorigin, _torigin);

  // Allocate memory for tempory data
  data = new char[_x*_y*_z*_t*_bytes];

  // Get pointer to image data
  ptr1 = _output->GetPointerToVoxels();

  // Convert data
  switch (_type) {
  case VOXEL_CHAR: {
    char *ptr2 = &(data[0]);

    // Read data
    this->ReadAsChar((char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<short>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_CHAR: {
    unsigned char *ptr2 = (unsigned char *)&(data[0]);

    // Read data
    this->ReadAsUChar((unsigned char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<short>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data
    this->ReadAsShort((short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<short>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_SHORT: {
    unsigned short *ptr2 = (unsigned short *)&(data[0]);

    // Read data
    this->ReadAsUShort((unsigned short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<short>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_INT: {
    int *ptr2 = (int *)&(data[0]);

    // Read data
    this->ReadAsInt((int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<short>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_INT: {
    unsigned int *ptr2 = (unsigned int *)&(data[0]);

    // Read data
    this->ReadAsUInt((unsigned int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<short>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data
    this->ReadAsFloat((float *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<short>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data
    this->ReadAsDouble((double *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<short>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_CHAR: {
    char *ptr2 = (char *)&(data[0]);

    // Read data.
    this->ReadAsChar((char *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<short>(*ptr2);
            ptr1->_y = static_cast<short>(*(ptr2 + 1));
            ptr1->_z = static_cast<short>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data.
    this->ReadAsShort((short *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<short>(*ptr2);
            ptr1->_y = static_cast<short>(*(ptr2 + 1));
            ptr1->_z = static_cast<short>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data.
    this->ReadAsFloat((float *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<short>(*ptr2);
            ptr1->_y = static_cast<short>(*(ptr2 + 1));
            ptr1->_z = static_cast<short>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data.
    this->ReadAsDouble((double *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<short>(*ptr2);
            ptr1->_y = static_cast<short>(*(ptr2 + 1));
            ptr1->_z = static_cast<short>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  default:
    cout << "unknown" << endl;
  }

  // Delete data
  delete []data;

  // Reflect if necessary
  if (_reflectX == True) _output->ReflectX();
  if (_reflectY == True) _output->ReflectY();
  if (_reflectZ == True) _output->ReflectZ();
}

template <> void irtkFileToImage<irtkVector3D<float> >::Run()
{
  char *data;
  int x, y, z, t;
  irtkVector3D<float> *ptr1;

  // Bring image to correct size
  _output->Initialize(_x, _y, _z, _t, _xsize, _ysize, _zsize, _tsize, _xaxis, _yaxis, _zaxis);
  _output->PutOrigin(_xorigin, _yorigin, _zorigin, _torigin);

  // Allocate memory for tempory data
  data = new char[_x*_y*_z*_t*_bytes];

  // Get pointer to image data
  ptr1 = _output->GetPointerToVoxels();

  // Convert data
  switch (_type) {
  case VOXEL_CHAR: {
    char *ptr2 = &(data[0]);

    // Read data
    this->ReadAsChar((char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<float>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_CHAR: {
    unsigned char *ptr2 = (unsigned char *)&(data[0]);

    // Read data
    this->ReadAsUChar((unsigned char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<float>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data
    this->ReadAsShort((short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<float>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_SHORT: {
    unsigned short *ptr2 = (unsigned short *)&(data[0]);

    // Read data
    this->ReadAsUShort((unsigned short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<float>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_INT: {
    int *ptr2 = (int *)&(data[0]);

    // Read data
    this->ReadAsInt((int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<float>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_INT: {
    unsigned int *ptr2 = (unsigned int *)&(data[0]);

    // Read data
    this->ReadAsUInt((unsigned int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<float>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data
    this->ReadAsFloat((float *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<float>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data
    this->ReadAsDouble((double *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<float>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_CHAR: {
    char *ptr2 = (char *)&(data[0]);

    // Read data.
    this->ReadAsChar((char *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<float>(*ptr2);
            ptr1->_y = static_cast<float>(*(ptr2 + 1));
            ptr1->_z = static_cast<float>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data.
    this->ReadAsShort((short *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<float>(*ptr2);
            ptr1->_y = static_cast<float>(*(ptr2 + 1));
            ptr1->_z = static_cast<float>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data.
    this->ReadAsFloat((float *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<float>(*ptr2);
            ptr1->_y = static_cast<float>(*(ptr2 + 1));
            ptr1->_z = static_cast<float>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data.
    this->ReadAsDouble((double *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<float>(*ptr2);
            ptr1->_y = static_cast<float>(*(ptr2 + 1));
            ptr1->_z = static_cast<float>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  default:
    cout << "unknown" << endl;
  }

  // Delete data
  delete []data;

  // Reflect if necessary
  if (_reflectX == True) _output->ReflectX();
  if (_reflectY == True) _output->ReflectY();
  if (_reflectZ == True) _output->ReflectZ();
}

template <> void irtkFileToImage<irtkVector3D<double> >::Run()
{
  char *data;
  int x, y, z, t;
  irtkVector3D<double> *ptr1;

  // Bring image to correct size
  _output->Initialize(_x, _y, _z, _t, _xsize, _ysize, _zsize, _tsize, _xaxis, _yaxis, _zaxis);
  _output->PutOrigin(_xorigin, _yorigin, _zorigin, _torigin);

  // Allocate memory for tempory data
  data = new char[_x*_y*_z*_t*_bytes];

  // Get pointer to image data
  ptr1 = _output->GetPointerToVoxels();

  // Convert data
  switch (_type) {
  case VOXEL_CHAR: {
    char *ptr2 = &(data[0]);

    // Read data
    this->ReadAsChar((char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<double>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_CHAR: {
    unsigned char *ptr2 = (unsigned char *)&(data[0]);

    // Read data
    this->ReadAsUChar((unsigned char *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<double>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data
    this->ReadAsShort((short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<double>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_SHORT: {
    unsigned short *ptr2 = (unsigned short *)&(data[0]);

    // Read data
    this->ReadAsUShort((unsigned short *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<double>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_INT: {
    int *ptr2 = (int *)&(data[0]);

    // Read data
    this->ReadAsInt((int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<double>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_U_INT: {
    unsigned int *ptr2 = (unsigned int *)&(data[0]);

    // Read data
    this->ReadAsUInt((unsigned int *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<double>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data
    this->ReadAsFloat((float *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<double>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data
    this->ReadAsDouble((double *)data, _x*_y*_z*_t, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = ptr1->_y = ptr1->_z = static_cast<double>(*ptr2);
            ptr1++;
            ptr2++;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_CHAR: {
    char *ptr2 = (char *)&(data[0]);

    // Read data.
    this->ReadAsChar((char *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<double>(*ptr2);
            ptr1->_y = static_cast<double>(*(ptr2 + 1));
            ptr1->_z = static_cast<double>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_SHORT: {
    short *ptr2 = (short *)&(data[0]);

    // Read data.
    this->ReadAsShort((short *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<double>(*ptr2);
            ptr1->_y = static_cast<double>(*(ptr2 + 1));
            ptr1->_z = static_cast<double>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_FLOAT: {
    float *ptr2 = (float *)&(data[0]);

    // Read data.
    this->ReadAsFloat((float *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<double>(*ptr2);
            ptr1->_y = static_cast<double>(*(ptr2 + 1));
            ptr1->_z = static_cast<double>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  case VOXEL_VECTOR_3D_DOUBLE: {
    double *ptr2 = (double *)&(data[0]);

    // Read data.
    this->ReadAsDouble((double *)data, _x*_y*_z*_t*3, _addr[0]);

    // Copy data
    for (t = 0; t < _t; t++) {
      for (z = 0; z < _z; z++) {
        for (y = 0; y < _y; y++) {
          for (x = 0; x < _x; x++) {
            ptr1->_x = static_cast<double>(*ptr2);
            ptr1->_y = static_cast<double>(*(ptr2 + 1));
            ptr1->_z = static_cast<double>(*(ptr2 + 2));

            ptr1++;
            ptr2 += 3;
          }
        }
      }
    }
  }
  break;
  default:
    cout << "unknown" << endl;
  }

  // Delete data
  delete []data;

  // Reflect if necessary
  if (_reflectX == True) _output->ReflectX();
  if (_reflectY == True) _output->ReflectY();
  if (_reflectZ == True) _output->ReflectZ();
}

template <class VoxelType> void irtkFileToImage<VoxelType>::Print()
{
  cout << "Name of class is " << this->NameOfClass() << endl;
  cout << "File name is " << _imagename << endl;
  cout << "Image dimensions are " << _x << " " << _y << " " << _z << " " << _t << endl;
  cout << "Image has " << _bytes << " bytes per voxel" << endl;
  cout << "Voxel dimensions are " << _xsize << " " << _ysize << " "
       << _zsize << " " << _tsize << endl;
  cout << "Voxel type is ";
  switch (_type) {
  case VOXEL_CHAR:
    cout << "char" << endl;
    break;
  case VOXEL_U_CHAR:
    cout << "unsigned char" << endl;
    break;
  case VOXEL_SHORT:
    cout << "short" << endl;
    break;
  case VOXEL_U_SHORT:
    cout << "unsigned short" << endl;
    break;
  case VOXEL_INT:
    cout << "int" << endl;
    break;
  case VOXEL_U_INT:
    cout << "unsigned int" << endl;
    break;
  case VOXEL_FLOAT:
    cout << "float" << endl;
    break;
  case VOXEL_DOUBLE:
    cout << "double" << endl;
    break;
  case VOXEL_VECTOR_3D_CHAR:
    cout << "vector 3d char" << endl;
    break;
  case VOXEL_VECTOR_3D_SHORT:
    cout << "vector 3d short" << endl;
    break;
  case VOXEL_VECTOR_3D_FLOAT:
    cout << "vector 3d float" << endl;
    break;
  case VOXEL_VECTOR_3D_DOUBLE:
    cout << "vector 3d double" << endl;
    break;
  default:
    cout << "unknown" << endl;
  }
}

template <class VoxelType> void irtkFileToImage<VoxelType>::Debug(char *message)
{
  if (_debug) cerr << message << endl;
}

template class irtkFileToImage<irtkBytePixel>;
template class irtkFileToImage<irtkGreyPixel>;
template class irtkFileToImage<irtkRealPixel>;
template class irtkFileToImage<irtkVector3D<char> >;
template class irtkFileToImage<irtkVector3D<short> >;
template class irtkFileToImage<irtkVector3D<float> >;
template class irtkFileToImage<irtkVector3D<double> >;
