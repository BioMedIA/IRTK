/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#define _IMPLEMENTS_GENERICIMAGE_

#include <irtkImage.h>

#include <irtkFileToImage.h>
#include <irtkImageToFile.h>

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, int t)
{
  _matrix = NULL;
  this->Initialize(x, y, z, t);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, double dx, double dy, double dz)
{
  _matrix = NULL;
  this->Initialize(x, y, z, dx, dy, dz);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, int t, double dx, double dy, double dz, double dt)
{
  _matrix = NULL;
  this->Initialize(x, y, z, t, dx, dy, dz, dt);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, double dx, double dy, double dz, const double *xaxis, const double *yaxis, const double *zaxis)
{
  _matrix = NULL;
  this->Initialize(x, y, z, dx, dy, dz, xaxis, yaxis, zaxis);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, int t, double dx, double dy, double dz, double dt, const double *xaxis, const double *yaxis, const double *zaxis)
{
  _matrix = NULL;
  this->Initialize(x, y, z, t, dx, dy, dz, dt, xaxis, yaxis, zaxis);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, double dx, double dy, double dz, irtkPoint origin, const double *xaxis, const double *yaxis, const double *zaxis)
{
  _matrix = NULL;
  this->Initialize(x, y, z, dx, dy, dz, origin, xaxis, yaxis, zaxis);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, int t, double dx, double dy, double dz, double dt, irtkPoint origin, double torigin, const double *xaxis, const double *yaxis, const double *zaxis)
{
  _matrix = NULL;
  this->Initialize(x, y, z, t, dx, dy, dz, dt, origin, torigin, xaxis, yaxis, zaxis);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(char *filename)
{
  _matrix = NULL;
  // Read image
  Read(filename);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(void)
{
  _matrix = NULL;
}

template <class VoxelType> irtkGenericImage<VoxelType>::~irtkGenericImage(void)
{
  if (_matrix != NULL) {
    Deallocate<VoxelType>(_matrix);
    _matrix = NULL;
  }
}

template <> const char *irtkGenericImage<irtkBytePixel>::NameOfClass()
{
  return "irtkGenericImage<BytePixel>";
}

template <> const char *irtkGenericImage<irtkGreyPixel>::NameOfClass()
{
  return "irtkGenericImage<GreyPixel>";
}

template <> const char *irtkGenericImage<irtkRealPixel>::NameOfClass()
{
  return "irtkGenericImage<RealPixel>";
}

template <> const char *irtkGenericImage<irtkRGBPixel>::NameOfClass()
{
  return "irtkGenericImage<RGBPixel>";
}

template <> const char *irtkGenericImage<irtkVector3D<char> >::NameOfClass()
{
  return "irtkGenericImage<irtkVector3D<char> >";
}

template <> const char *irtkGenericImage<irtkVector3D<short> >::NameOfClass()
{
  return "irtkGenericImage<irtkVector3D<short> >";
}

template <> const char *irtkGenericImage<irtkVector3D<float> >::NameOfClass()
{
  return "irtkGenericImage<irtkVector3D<float> >";
}

template <> const char *irtkGenericImage<irtkVector3D<double> >::NameOfClass()
{
  return "irtkGenericImage<irtkVector3D<double> >";
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize(int x, int y, int z, int t)
{
  int i, n;
  VoxelType *ptr;

  // Free memory
  if ((_x != x) || (_y != y) || (_z != z) || (_t != t)) {
    // Free old memory
    if (_matrix != NULL) Deallocate<VoxelType>(_matrix);
    // Allocate new memory
    _matrix = Allocate(_matrix, x, y, z, t);
  }

  // Initialize base class
  this->irtkBaseImage::Initialize(x, y, z, t);

  // Initialize voxels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType();
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize(int x, int y, int z, double dx, double dy, double dz)
{
  int i, n;
  VoxelType *ptr;

  // Free memory
  if ((_x != x) || (_y != y) || (_z != z) || (_t != 1)) {
    // Free old memory
    if (_matrix != NULL) Deallocate<VoxelType>(_matrix);
    // Allocate new memory
    _matrix = Allocate(_matrix, x, y, z, 1);
  }

  // Initialize base class
  this->irtkBaseImage::Initialize(x, y, z, dx, dy, dz);

  // Initialize voxels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType();
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize(int x, int y, int z, int t, double dx, double dy, double dz, double dt)
{
  int i, n;
  VoxelType *ptr;

  // Free memory
  if ((_x != x) || (_y != y) || (_z != z) || (_t != t)) {
    // Free old memory
    if (_matrix != NULL) Deallocate<VoxelType>(_matrix);
    // Allocate new memory
    _matrix = Allocate(_matrix, x, y, z, t);
  }

  // Initialize base class
  this->irtkBaseImage::Initialize(x, y, z, t, dx, dy, dz, dt);

  // Initialize voxels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType();
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize(int x, int y, int z, double dx, double dy, double dz, const double *xaxis, const double *yaxis, const double *zaxis)
{
  int i, n;
  VoxelType *ptr;

  // Free memory
  if ((_x != x) || (_y != y) || (_z != z) || (_t != 1)) {
    // Free old memory
    if (_matrix != NULL) Deallocate<VoxelType>(_matrix);
    // Allocate new memory
    _matrix = Allocate(_matrix, x, y, z, 1);
  }

  // Initialize base class
  this->irtkBaseImage::Initialize(x, y, z, dx, dy, dz, xaxis, yaxis, zaxis);

  // Initialize voxels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType();
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize(int x, int y, int z, int t, double dx, double dy, double dz, double dt, const double *xaxis, const double *yaxis, const double *zaxis)
{
  int i, n;
  VoxelType *ptr;

  // Free memory
  if ((_x != x) || (_y != y) || (_z != z) || (_t != t)) {
    // Free old memory
    if (_matrix != NULL) Deallocate<VoxelType>(_matrix);
    // Allocate new memory
    _matrix = Allocate(_matrix, x, y, z, t);
  }

  // Initialize base class
  this->irtkBaseImage::Initialize(x, y, z, t, dx, dy, dz, dt, xaxis, yaxis, zaxis);

  // Initialize voxels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType();
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize(int x, int y, int z, double dx, double dy, double dz, irtkPoint origin, const double *xaxis, const double *yaxis, const double *zaxis)
{
  int i, n;
  VoxelType *ptr;

  // Free memory
  if ((_x != x) || (_y != y) || (_z != z) || (_t != 1)) {
    // Free old memory
    if (_matrix != NULL) Deallocate<VoxelType>(_matrix);
    // Allocate new memory
    _matrix = Allocate(_matrix, x, y, z, 1);
  }

  // Initialize base class
  this->irtkBaseImage::Initialize(x, y, z, dx, dy, dz, origin, xaxis, yaxis, zaxis);

  // Initialize voxels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType();
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize(int x, int y, int z, int t, double dx, double dy, double dz, double dt, irtkPoint origin, double torigin, const double *xaxis, const double *yaxis, const double *zaxis)
{
  int i, n;
  VoxelType *ptr;

  // Free memory
  if ((_x != x) || (_y != y) || (_z != z) || (_t != t)) {
    // Free old memory
    if (_matrix != NULL) Deallocate<VoxelType>(_matrix);
    // Allocate new memory
    _matrix = Allocate(_matrix, x, y, z, t);
  }

  // Initialize base class
  this->irtkBaseImage::Initialize(x, y, z, t, dx, dy, dz, dt, origin, torigin, xaxis, yaxis, zaxis);

  // Initialize voxels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType();
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize(const irtkBaseImage &image)
{
  int i, n;
  VoxelType *ptr;

  // Free memory
  if ((_x != image.GetX()) || (_y != image.GetY()) || (_z != image.GetZ()) || (_t != image.GetT())) {
    // Free old memory
    if (_matrix != NULL) Deallocate<VoxelType>(_matrix);
    // Allocate new memory
    _matrix = Allocate(_matrix, image.GetX(), image.GetY(), image.GetZ(), image.GetT());
  }

  // Initialize base class
  this->irtkBaseImage::Initialize(image);

  // Initialize voxels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType();
  }
}

template <> void irtkGenericImage<irtkRGBPixel>::Read(const char *filename)
{}

template <class VoxelType> void irtkGenericImage<VoxelType>::Read(const char *filename)
{
  // Allocate file reader
  irtkFileToImage<VoxelType> *reader =
    irtkFileToImage<VoxelType>::New(filename);

  // Set output
  reader->SetOutput(this);

  // Run reader
  reader->Run();

  // Delete reader
  delete reader;
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Write(const char *filename)
{
  // Allocate file reader
  irtkImageToFile<VoxelType> *writer =
    irtkImageToFile<VoxelType>::New(filename);

  // Set input
  writer->SetInput(this);

  // Run writer
  writer->Run();

  // Delete writer
  delete writer;
}

template <class VoxelType> void irtkGenericImage<VoxelType>::GetMinMax(VoxelType *min, VoxelType *max) const
{
  int i, n;
  VoxelType *ptr;

  *min = VoxelType();
  *max = VoxelType();

  // Initialize pixels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();

  if (n > 0) {
    *min = ptr[0];
    *max = ptr[0];
    for (i = 0; i < n; i++) {
      if (ptr[i] < *min) *min = ptr[i];
      if (ptr[i] > *max) *max = ptr[i];
    }
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::GetMinMaxPad(VoxelType *min, VoxelType *max, VoxelType pad) const
{
  int i, n;
  VoxelType *ptr;
  bool first=true;

  *min = VoxelType();
  *max = VoxelType();

  // Initialize pixels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();

  if (n > 0) {
    for (i = 0; i < n; i++) {
      if (ptr[i]!=pad) {
        if (first) {
          first=false;
          *min = ptr[i];
          *max = ptr[i];
        } else {
          if (ptr[i] < *min) *min = ptr[i];
          if (ptr[i] > *max) *max = ptr[i];
        }
      }
    }
  }
}

template <> void irtkGenericImage<irtkBytePixel>::PutMinMax(irtkBytePixel min, irtkBytePixel max)
{
  int i, n;
  irtkBytePixel *ptr, min_val, max_val;

  // Get lower and upper bound
  this->GetMinMax(&min_val, &max_val);

  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = round(((ptr[i] - min_val) / double(max_val - min_val)) *
                   (max - min) + min);
  }
}

template <> void irtkGenericImage<irtkGreyPixel>::PutMinMax(irtkGreyPixel min, irtkGreyPixel max)
{
  int i, n;
  irtkGreyPixel *ptr, min_val, max_val;

  // Get lower and upper bound
  this->GetMinMax(&min_val, &max_val);

  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = round(((ptr[i] - min_val) / double(max_val - min_val)) *
                   (max - min) + min);
  }
}

template <> void irtkGenericImage<irtkVector3D<char> >::PutMinMax(irtkVector3D<char> min, irtkVector3D<char> max)
{
  irtkVector3D<char> min_val, max_val;
  min_val._x = CHAR_MAX;
  min_val._y = CHAR_MAX;
  min_val._z = CHAR_MAX;
  max_val._x = CHAR_MIN;
  max_val._y = CHAR_MIN;
  max_val._z = CHAR_MIN;

  // Get lower and upper bound of components.
  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          if (_matrix[t][z][y][x]._x < min_val._x) {
            min_val._x = _matrix[t][z][y][x]._x;
          }
          if (_matrix[t][z][y][x]._y < min_val._y) {
            min_val._y = _matrix[t][z][y][x]._y;
          }
          if (_matrix[t][z][y][x]._z < min_val._z) {
            min_val._z = _matrix[t][z][y][x]._z;
          }
          if (_matrix[t][z][y][x]._x > max_val._x) {
            max_val._x = _matrix[t][z][y][x]._x;
          }
          if (_matrix[t][z][y][x]._y > max_val._y) {
            max_val._y = _matrix[t][z][y][x]._y;
          }
          if (_matrix[t][z][y][x]._z > max_val._z) {
            max_val._z = _matrix[t][z][y][x]._z;
          }
        }
      }
    }
  }

  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          _matrix[t][z][y][x]._x = static_cast<char>((_matrix[t][z][y][x]._x - min_val._x)*(max._x - min._x)/static_cast<double>(max_val._x - min_val._x) + min._x);
          _matrix[t][z][y][x]._y = static_cast<char>((_matrix[t][z][y][x]._y - min_val._y)*(max._y - min._y)/static_cast<double>(max_val._y - min_val._y) + min._y);
          _matrix[t][z][y][x]._z = static_cast<char>((_matrix[t][z][y][x]._z - min_val._z)*(max._z - min._z)/static_cast<double>(max_val._z - min_val._z) + min._z);
        }
      }
    }
  }
}

template <> void irtkGenericImage<irtkVector3D<short> >::PutMinMax(irtkVector3D<short> min, irtkVector3D<short> max)
{
  irtkVector3D<short> min_val, max_val;
  min_val._x = SHRT_MAX;
  min_val._y = SHRT_MAX;
  min_val._z = SHRT_MAX;
  max_val._x = SHRT_MIN;
  max_val._y = SHRT_MIN;
  max_val._z = SHRT_MIN;

  // Get lower and upper bound of components.
  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          if (_matrix[t][z][y][x]._x < min_val._x) {
            min_val._x = _matrix[t][z][y][x]._x;
          }
          if (_matrix[t][z][y][x]._y < min_val._y) {
            min_val._y = _matrix[t][z][y][x]._y;
          }
          if (_matrix[t][z][y][x]._z < min_val._z) {
            min_val._z = _matrix[t][z][y][x]._z;
          }
          if (_matrix[t][z][y][x]._x > max_val._x) {
            max_val._x = _matrix[t][z][y][x]._x;
          }
          if (_matrix[t][z][y][x]._y > max_val._y) {
            max_val._y = _matrix[t][z][y][x]._y;
          }
          if (_matrix[t][z][y][x]._z > max_val._z) {
            max_val._z = _matrix[t][z][y][x]._z;
          }
        }
      }
    }
  }

  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          _matrix[t][z][y][x]._x = static_cast<short>((_matrix[t][z][y][x]._x - min_val._x)*(max._x - min._x)/static_cast<double>(max_val._x - min_val._x) + min._x);
          _matrix[t][z][y][x]._y = static_cast<short>((_matrix[t][z][y][x]._y - min_val._y)*(max._y - min._y)/static_cast<double>(max_val._y - min_val._y) + min._y);
          _matrix[t][z][y][x]._z = static_cast<short>((_matrix[t][z][y][x]._z - min_val._z)*(max._z - min._z)/static_cast<double>(max_val._z - min_val._z) + min._z);
        }
      }
    }
  }
}

template <> void irtkGenericImage<irtkVector3D<float> >::PutMinMax(irtkVector3D<float> min, irtkVector3D<float> max)
{
  irtkVector3D<float> min_val, max_val;
  min_val._x = FLT_MAX;
  min_val._y = FLT_MAX;
  min_val._z = FLT_MAX;
  max_val._x = FLT_MIN;
  max_val._y = FLT_MIN;
  max_val._z = FLT_MIN;

  // Get lower and upper bound of components.
  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          if (_matrix[t][z][y][x]._x < min_val._x) {
            min_val._x = _matrix[t][z][y][x]._x;
          }
          if (_matrix[t][z][y][x]._y < min_val._y) {
            min_val._y = _matrix[t][z][y][x]._y;
          }
          if (_matrix[t][z][y][x]._z < min_val._z) {
            min_val._z = _matrix[t][z][y][x]._z;
          }
          if (_matrix[t][z][y][x]._x > max_val._x) {
            max_val._x = _matrix[t][z][y][x]._x;
          }
          if (_matrix[t][z][y][x]._y > max_val._y) {
            max_val._y = _matrix[t][z][y][x]._y;
          }
          if (_matrix[t][z][y][x]._z > max_val._z) {
            max_val._z = _matrix[t][z][y][x]._z;
          }
        }
      }
    }
  }

  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          _matrix[t][z][y][x]._x = static_cast<float>((_matrix[t][z][y][x]._x - min_val._x)*(max._x - min._x)/static_cast<double>(max_val._x - min_val._x) + min._x);
          _matrix[t][z][y][x]._y = static_cast<float>((_matrix[t][z][y][x]._y - min_val._y)*(max._y - min._y)/static_cast<double>(max_val._y - min_val._y) + min._y);
          _matrix[t][z][y][x]._z = static_cast<float>((_matrix[t][z][y][x]._z - min_val._z)*(max._z - min._z)/static_cast<double>(max_val._z - min_val._z) + min._z);
        }
      }
    }
  }
}

template <> void irtkGenericImage<irtkVector3D<double> >::PutMinMax(irtkVector3D<double> min, irtkVector3D<double> max)
{
  irtkVector3D<double> min_val, max_val;
  min_val._x = DBL_MAX;
  min_val._y = DBL_MAX;
  min_val._z = DBL_MAX;
  max_val._x = DBL_MIN;
  max_val._y = DBL_MIN;
  max_val._z = DBL_MIN;

  // Get lower and upper bound of components.
  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          if (_matrix[t][z][y][x]._x < min_val._x) {
            min_val._x = _matrix[t][z][y][x]._x;
          }
          if (_matrix[t][z][y][x]._y < min_val._y) {
            min_val._y = _matrix[t][z][y][x]._y;
          }
          if (_matrix[t][z][y][x]._z < min_val._z) {
            min_val._z = _matrix[t][z][y][x]._z;
          }
          if (_matrix[t][z][y][x]._x > max_val._x) {
            max_val._x = _matrix[t][z][y][x]._x;
          }
          if (_matrix[t][z][y][x]._y > max_val._y) {
            max_val._y = _matrix[t][z][y][x]._y;
          }
          if (_matrix[t][z][y][x]._z > max_val._z) {
            max_val._z = _matrix[t][z][y][x]._z;
          }
        }
      }
    }
  }

  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          _matrix[t][z][y][x]._x = static_cast<double>((_matrix[t][z][y][x]._x - min_val._x)*(max._x - min._x)/static_cast<double>(max_val._x - min_val._x) + min._x);
          _matrix[t][z][y][x]._y = static_cast<double>((_matrix[t][z][y][x]._y - min_val._y)*(max._y - min._y)/static_cast<double>(max_val._y - min_val._y) + min._y);
          _matrix[t][z][y][x]._z = static_cast<double>((_matrix[t][z][y][x]._z - min_val._z)*(max._z - min._z)/static_cast<double>(max_val._z - min_val._z) + min._z);
        }
      }
    }
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::PutMinMax(VoxelType min, VoxelType max)
{
  int i, n;
  VoxelType *ptr, min_val, max_val;

  // Get lower and upper bound
  this->GetMinMax(&min_val, &max_val);

  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType(((ptr[i] - min_val) / (max_val - min_val)) *
                       (max - min) + min);
  }
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::GetRegion(int k, int m) const
{
  int i, j;
  double x1, y1, z1, t1, x2, y2, z2, t2;

  if ((k < 0) || (k >= _z) || (m < 0) || (m >= _t)) {
    cerr << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
    exit(1);
  }

  // Init
  irtkGenericImage<VoxelType> image(_x, _y, 1, 1, _dx, _dy, _dz, _dt, _xaxis, _yaxis, _zaxis);

  // Calculate position of first voxel in roi in original image
  x1 = 0;
  y1 = 0;
  z1 = k;
  this->ImageToWorld(x1, y1, z1);
  t1 = this->ImageToTime(m);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  t2 = 0;
  image.ImageToWorld(x2, y2, z2);
  t2 = image.ImageToTime(0);

  // Shift origin of new image accordingly
  image.PutOrigin(x1 - x2, y1 - y2, z1 - z2, t1 - t2);

  // Copy region
  for (j = 0; j < _y; j++) {
    for (i = 0; i < _x; i++) {
      image._matrix[0][0][j][i] = _matrix[m][k][j][i];
    }
  }
  return image;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::GetFrame(int l) const
{
  int i, j, k;

  if ((l < 0) || (l >= _t)) {
    cerr << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
    exit(1);
  }

  // Init
  irtkGenericImage<VoxelType> image(_x, _y, _z, 1, _dx, _dy, _dz, _dt, _origin, _torigin, _xaxis, _yaxis, _zaxis);
  image.PutOrigin(_origin);

  // Copy region
  for (k = 0; k < _z; k++) {
    for (j = 0; j < _y; j++) {
      for (i = 0; i < _x; i++) {
        image._matrix[0][k][j][i] = _matrix[l][k][j][i];
      }
    }
  }
  return image;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::GetRegion(int i1, int j1, int k1, int i2, int j2, int k2) const
{
  int i, j, k, l;
  double x1, y1, z1, x2, y2, z2;

  if ((i1 < 0) || (i1 >= i2) ||
      (j1 < 0) || (j1 >= j2) ||
      (k1 < 0) || (k1 >= k2) ||
      (i2 > _x) || (j2 > _y) || (k2 > _z)) {
    cerr << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
    exit(1);
  }

  // Init
  irtkGenericImage<VoxelType> image(i2 - i1, j2 - j1, k2 - k1, _t, _dx, _dy, _dz, _dt, _xaxis, _yaxis, _zaxis);

  // Calculate position of first voxel in roi in original image
  x1 = i1;
  y1 = j1;
  z1 = k1;
  this->ImageToWorld(x1, y1, z1);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  image.ImageToWorld(x2, y2, z2);

  // Shift origin of new image accordingly
  image.PutOrigin(x1 - x2, y1 - y2, z1 - z2);

  // Copy region
  for (l = 0; l < _t; l++) {
    for (k = k1; k < k2; k++) {
      for (j = j1; j < j2; j++) {
        for (i = i1; i < i2; i++) {
          image._matrix[l][k-k1][j-j1][i-i1] = _matrix[l][k][j][i];
        }
      }
    }
  }
  return image;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::GetRegion(int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2) const
{
  int i, j, k, l;
  double x1, y1, z1, x2, y2, z2;

  if ((i1 < 0) || (i1 >= i2) ||
      (j1 < 0) || (j1 >= j2) ||
      (k1 < 0) || (k1 >= k2) ||
      (l1 < 0) || (l1 >= l2) ||
      (i2 > _x) || (j2 > _y) || (k2 > _z) || (l2 > _t)) {
    cerr << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
    exit(1);
  }

  // Init
  irtkGenericImage<VoxelType> image(i2 - i1, j2 - j1, k2 - k1, l2 - l1, _dx, _dy, _dz, _dt, _xaxis, _yaxis, _zaxis);

  // Calculate position of first voxel in roi in original image
  x1 = i1;
  y1 = j1;
  z1 = k1;
  this->ImageToWorld(x1, y1, z1);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  image.ImageToWorld(x2, y2, z2);

  // Shift origin of new image accordingly
  image.PutOrigin(x1 - x2, y1 - y2, z1 - z2);

  // Copy region
  for (l = l1; l < l2; l++) {
    for (k = k1; k < k2; k++) {
      for (j = j1; j < j2; j++) {
        for (i = i1; i < i2; i++) {
          image._matrix[l-l1][k-k1][j-j1][i-i1] = _matrix[l][k][j][i];
        }
      }
    }
  }
  return image;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator+=(const irtkGenericImage<VoxelType> &image)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  if (! (this->irtkBaseImage::operator==(image))) {
    cerr << "irtkGenericImage<VoxelType>::operator+=: Size mismatch in images\n";
    exit(1);
  }

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] += ptr2[i];
  }
  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator+(const irtkGenericImage<VoxelType> &image)
{
  irtkGenericImage<VoxelType> tmp(*this); tmp += image; return tmp;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator-=(const irtkGenericImage<VoxelType> &image)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  if (! (this->irtkBaseImage::operator==(image))) {
    cerr << "irtkGenericImage<VoxelType>::operator-=: Size mismatch in images\n";
    exit(1);
  }

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] -= ptr2[i];
  }
  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator-(const irtkGenericImage<VoxelType> &image)
{
  irtkGenericImage<VoxelType> tmp(*this); tmp -= image; return tmp;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator*=(const irtkGenericImage<VoxelType> &image)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  if (! (this->irtkBaseImage::operator==(image))) {
    cerr << "irtkGenericImage<VoxelType>::operator*=: Size mismatch in images\n";
    exit(1);
  }

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] *= ptr2[i];
  }
  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator*(const irtkGenericImage<VoxelType> &image)
{
  irtkGenericImage<VoxelType> tmp(*this); tmp *= image; return tmp;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator/=(const irtkGenericImage<VoxelType> &image)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  if (! (this->irtkBaseImage::operator==(image))) {
    cerr << "irtkGenericImage<VoxelType>::operator/=: Size mismatch in images\n";
    exit(1);
  }

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if (ptr2[i] != VoxelType()) {
      ptr1[i] /= ptr2[i];
    } else {
      ptr1[i] = VoxelType();
    }
  }
  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator/(const irtkGenericImage<VoxelType> &image)
{
  irtkGenericImage<VoxelType> tmp(*this); tmp /= image; return tmp;
}

template <class VoxelType> Bool irtkGenericImage<VoxelType>::operator==(const irtkGenericImage<VoxelType> &image)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  if (! (this->irtkBaseImage::operator==(image))) {
    return False;
  }

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if (ptr1[i] != ptr2[i]) return False;
  }
  return True;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator+=(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr;

  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] += pixel;
  }
  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator+(VoxelType pixel)
{
  irtkGenericImage<VoxelType> tmp(*this); tmp += pixel; return tmp;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator-=(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr;

  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] -= pixel;
  }
  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator-(VoxelType pixel)
{
  irtkGenericImage<VoxelType> tmp(*this); tmp -= pixel; return tmp;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator*=(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr;

  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] *= pixel;
  }
  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator*(VoxelType pixel)
{
  irtkGenericImage<VoxelType> tmp(*this); tmp *= pixel; return tmp;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator/=(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr;

  if (pixel != VoxelType()) {
    n   = this->GetNumberOfVoxels();
    ptr = this->GetPointerToVoxels();
    for (i = 0; i < n; i++) {
      ptr[i] /= pixel;
    }
  } else {
    cerr << "irtkGenericImage<VoxelType>::operator/=: Division by zero" << endl;
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<char> >& irtkGenericImage<irtkVector3D<char> >::operator/=(irtkVector3D<char> pixel)
{
  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          _matrix[t][z][y][x]._x /= pixel._x;
          _matrix[t][z][y][x]._y /= pixel._y;
          _matrix[t][z][y][x]._z /= pixel._z;
        }
      }
    }
  }

  return *this;
}

template <> irtkGenericImage<irtkVector3D<short> >& irtkGenericImage<irtkVector3D<short> >::operator/=(irtkVector3D<short> pixel)
{
  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          _matrix[t][z][y][x]._x /= pixel._x;
          _matrix[t][z][y][x]._y /= pixel._y;
          _matrix[t][z][y][x]._z /= pixel._z;
        }
      }
    }
  }

  return *this;
}

template <> irtkGenericImage<irtkVector3D<float> >& irtkGenericImage<irtkVector3D<float> >::operator/=(irtkVector3D<float> pixel)
{
  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          _matrix[t][z][y][x]._x /= pixel._x;
          _matrix[t][z][y][x]._y /= pixel._y;
          _matrix[t][z][y][x]._z /= pixel._z;
        }
      }
    }
  }

  return *this;
}

template <> irtkGenericImage<irtkVector3D<double> >& irtkGenericImage<irtkVector3D<double> >::operator/=(irtkVector3D<double> pixel)
{
  for (int t = 0; t < _t; t++) {
    for (int z = 0; z < _z; z++) {
      for (int y = 0; y < _y; y++) {
        for (int x = 0; x < _x; x++) {
          _matrix[t][z][y][x]._x /= pixel._x;
          _matrix[t][z][y][x]._y /= pixel._y;
          _matrix[t][z][y][x]._z /= pixel._z;
        }
      }
    }
  }

  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator/(VoxelType pixel)
{
  irtkGenericImage<VoxelType> tmp(*this); tmp /= pixel; return tmp;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator>(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  irtkGenericImage<VoxelType> image(*this);

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if (ptr1[i] > pixel) {
      ptr2[i] = pixel;
    } else {
      ptr2[i] = ptr1[i];
    }
  }
  return image;
}

template <class VoxelType> irtkGenericImage<VoxelType> &irtkGenericImage<VoxelType>::operator>=(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr;

  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if (ptr[i] > pixel) {
      ptr[i] = pixel;
    }
  }
  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator<(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  irtkGenericImage<VoxelType> image(*this);

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if (ptr1[i] < pixel) {
      ptr2[i] = pixel;
    } else {
      ptr2[i] = ptr1[i];
    }
  }
  return image;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator<=(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr;

  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if (ptr[i] < pixel) {
      ptr[i] = pixel;
    }
  }
  return *this;
}

template <class VoxelType> void irtkGenericImage<VoxelType>::ReflectX()
{
  int x, y, z, t;

  for (t = 0; t < _t; t++) {
    for (z = 0; z < _z; z++) {
      for (y = 0; y < _y; y++) {
        for (x = 0; x < _x / 2; x++) {
          swap(_matrix[t][z][y][x], _matrix[t][z][y][_x-(x+1)]);
        }
      }
    }
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::ReflectY()
{
  int x, y, z, t;

  for (t = 0; t < _t; t++) {
    for (z = 0; z < _z; z++) {
      for (y = 0; y < _y / 2; y++) {
        for (x = 0; x < _x; x++) {
          swap(_matrix[t][z][y][x], _matrix[t][z][_y-(y+1)][x]);
        }
      }
    }
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::ReflectZ()
{
  int x, y, z, t;

  for (t = 0; t < _t; t++) {
    for (z = 0; z < _z / 2; z++) {
      for (y = 0; y < _y; y++) {
        for (x = 0; x < _x; x++) {
          swap(_matrix[t][z][y][x], _matrix[t][_z-(z+1)][y][x]);
        }
      }
    }
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipXY()
{
  int i, j, k, m;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _y, _x, _z, _t);

  for (m = 0; m < _t; m++) {
    for (k = 0; k < _z; k++) {
      for (j = 0; j < _y; j++) {
        for (i = 0; i < _x; i++) {
          matrix[m][k][i][j] = _matrix[m][k][j][i];
        }
      }
    }
  }

  // Swap pointers
  swap(matrix, _matrix);

  // Deallocate memory
  matrix = Deallocate<VoxelType>(matrix);

  // Swap image dimensions
  swap(_x, _y);

  // Swap voxel dimensions
  swap(_dx, _dy);

  // Swap voxel dimensions
  swap(_origin._x, _origin._y);

  // Update transformation matrix
  this->UpdateMatrix();

}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipXZ()
{
  int i, j, k, l;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _z, _y, _x, _t);

  for (l = 0; l < _t; l++) {
    for (k = 0; k < _z; k++) {
      for (j = 0; j < _y; j++) {
        for (i = 0; i < _x; i++) {
          matrix[l][i][j][k] = _matrix[l][k][j][i];
        }
      }
    }
  }

  // Swap pointers
  swap(matrix, _matrix);

  // Deallocate memory
  matrix = Deallocate<VoxelType>(matrix);

  // Swap image dimensions
  swap(_x, _z);

  // Swap voxel dimensions
  swap(_dx, _dz);

  // Swap voxel dimensions
  swap(_origin._x, _origin._z);

  // Update transformation matrix
  this->UpdateMatrix();
}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipYZ()
{
  int i, j, k, l;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _x, _z, _y, _t);

  for (l = 0; l < _t; l++) {
    for (k = 0; k < _z; k++) {
      for (j = 0; j < _y; j++) {
        for (i = 0; i < _x; i++) {
          matrix[l][j][k][i] = _matrix[l][k][j][i];
        }
      }
    }
  }

  // Swap pointers
  swap(matrix, _matrix);

  // Deallocate memory
  matrix = Deallocate<VoxelType>(matrix);

  // Swap image dimensions
  swap(_y, _z);

  // Swap voxel dimensions
  swap(_dy, _dz);

  // Swap voxel dimensions
  swap(_origin._y, _origin._z);

  // Update transformation matrix
  this->UpdateMatrix();
}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipXT()
{
  int i, j, k, m;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _t, _y, _z, _x);

  for (m = 0; m < _t; m++) {
    for (k = 0; k < _z; k++) {
      for (j = 0; j < _y; j++) {
        for (i = 0; i < _x; i++) {
          matrix[i][k][j][m] = _matrix[m][k][j][i];
        }
      }
    }
  }

  // Swap pointers
  swap(matrix, _matrix);

  // Deallocate memory
  matrix = Deallocate<VoxelType>(matrix);

  // Swap image dimensions
  swap(_x, _t);

  // Swap voxel dimensions
  swap(_dx, _dt);

  // Swap voxel dimensions
  swap(_origin._x, _torigin);

  // Update transformation matrix
  this->UpdateMatrix();

}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipYT()
{
  int i, j, k, m;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _x, _t, _z, _y);

  for (m = 0; m < _t; m++) {
    for (k = 0; k < _z; k++) {
      for (j = 0; j < _y; j++) {
        for (i = 0; i < _x; i++) {
          matrix[j][k][m][i] = _matrix[m][k][j][i];
        }
      }
    }
  }

  // Swap pointers
  swap(matrix, _matrix);

  // Deallocate memory
  matrix = Deallocate<VoxelType>(matrix);

  // Swap image dimensions
  swap(_y, _t);

  // Swap voxel dimensions
  swap(_dy, _dt);

  // Swap voxel dimensions
  swap(_origin._y, _torigin);

  // Update transformation matrix
  this->UpdateMatrix();

}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipZT()
{
  int i, j, k, m;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _x, _y, _t, _z);

  for (m = 0; m < _t; m++) {
    for (k = 0; k < _z; k++) {
      for (j = 0; j < _y; j++) {
        for (i = 0; i < _x; i++) {
          matrix[k][m][j][i] = _matrix[m][k][j][i];
        }
      }
    }
  }

  // Swap pointers
  swap(matrix, _matrix);

  // Deallocate memory
  matrix = Deallocate<VoxelType>(matrix);

  // Swap image dimensions
  swap(_z, _t);

  // Swap voxel dimensions
  swap(_dz, _dt);

  // Swap voxel dimensions
  swap(_origin._z, _torigin);

  // Update transformation matrix
  this->UpdateMatrix();

}

template <> irtkGenericImage<irtkBytePixel>& irtkGenericImage<irtkBytePixel>::operator=
(const irtkGenericImage<irtkBytePixel> &image)
{
  int i, n;
  irtkBytePixel *ptr1, *ptr2;

  if (this == &image) return *this;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkBytePixel>& irtkGenericImage<irtkBytePixel>::operator=
(const irtkGenericImage<irtkGreyPixel> &image)
{
  int i, n;
  irtkBytePixel *ptr1;
  irtkGreyPixel *ptr2, min, max;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = round(MAX_BYTE*(ptr2[i] - min)/(double)(max - min));
  }
  return *this;
}

template <> irtkGenericImage<irtkBytePixel>& irtkGenericImage<irtkBytePixel>::operator=
(const irtkGenericImage<irtkRealPixel> &image)
{
  int i, n;
  irtkBytePixel *ptr1;
  irtkRealPixel *ptr2, min, max;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = round(MAX_BYTE*(ptr2[i] - min)/(double)(max - min));
  }
  return *this;
}

template <> irtkGenericImage<irtkBytePixel>& irtkGenericImage<irtkBytePixel>::operator=
(const irtkGenericImage<irtkVector3D<char> >& image)
{
  int i, n;
  irtkBytePixel* ptr1;
  irtkVector3D<char>* ptr2, min, max;
  double lmin, lmax;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round(MAX_BYTE*(l - lmin)/(lmax - lmin));
  }
  return *this;
}

template <> irtkGenericImage<irtkBytePixel>& irtkGenericImage<irtkBytePixel>::operator=
(const irtkGenericImage<irtkVector3D<short> >& image)
{
  int i, n;
  irtkBytePixel* ptr1;
  irtkVector3D<short>* ptr2, min, max;
  double lmin, lmax;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round(MAX_BYTE*(l - lmin)/(lmax - lmin));
  }
  return *this;
}

template <> irtkGenericImage<irtkBytePixel>& irtkGenericImage<irtkBytePixel>::operator=
(const irtkGenericImage<irtkVector3D<float> >& image)
{
  int i, n;
  irtkBytePixel* ptr1;
  irtkVector3D<float>* ptr2, min, max;
  double lmin, lmax;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round(MAX_BYTE*(l - lmin)/(lmax - lmin));
  }
  return *this;
}

template <> irtkGenericImage<irtkBytePixel>& irtkGenericImage<irtkBytePixel>::operator=
(const irtkGenericImage<irtkVector3D<double> >& image)
{
  int i, n;
  irtkBytePixel* ptr1;
  irtkVector3D<double>* ptr2, min, max;
  double lmin, lmax;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(min._x*min._x + min._y*min._y + min._z*min._z);
  lmax = sqrt(max._x*max._x + max._y*max._y + max._z*max._z);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round(MAX_BYTE*(l - lmin)/(lmax - lmin));
  }
  return *this;
}

template <>  irtkGenericImage<irtkGreyPixel>& irtkGenericImage<irtkGreyPixel>::operator=
(const irtkGenericImage<irtkBytePixel> &image)
{
  int i, n;
  irtkGreyPixel *ptr1;
  irtkBytePixel *ptr2;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkGreyPixel>& irtkGenericImage<irtkGreyPixel>::operator=
(const irtkGenericImage<irtkGreyPixel> &image)
{
  int i, n;
  irtkGreyPixel *ptr1, *ptr2;

  if (this == &image) return *this;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkGreyPixel>& irtkGenericImage<irtkGreyPixel>::operator=
(const irtkGenericImage<irtkRealPixel> &image)
{
  int i, n;
  irtkGreyPixel *ptr1;
  irtkRealPixel *ptr2, min, max;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = round(((ptr2[i] - min)/(double)(max - min)) *
                    (MAX_GREY - MIN_GREY) + MIN_GREY);
  }
  return *this;
}

template <> irtkGenericImage<irtkGreyPixel>& irtkGenericImage<irtkGreyPixel>::operator=
(const irtkGenericImage<irtkVector3D<char> >& image)
{
  int i, n;
  irtkGreyPixel* ptr1;
  irtkVector3D<char>* ptr2, min, max;
  double lmin, lmax;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round((l - lmin)/(double)(lmax - lmin)*(MAX_GREY - MIN_GREY) + MIN_GREY);
  }
  return *this;
}

template <> irtkGenericImage<irtkGreyPixel>& irtkGenericImage<irtkGreyPixel>::operator=
(const irtkGenericImage<irtkVector3D<short> >& image)
{
  int i, n;
  irtkGreyPixel* ptr1;
  irtkVector3D<short>* ptr2, min, max;
  double lmin, lmax;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round((l - lmin)/(double)(lmax - lmin)*(MAX_GREY - MIN_GREY) + MIN_GREY);
  }
  return *this;
}

template <> irtkGenericImage<irtkGreyPixel>& irtkGenericImage<irtkGreyPixel>::operator=
(const irtkGenericImage<irtkVector3D<float> >& image)
{
  int i, n;
  irtkGreyPixel* ptr1;
  irtkVector3D<float>* ptr2, min, max;
  double lmin, lmax;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round((l - lmin)/(double)(lmax - lmin)*(MAX_GREY - MIN_GREY) + MIN_GREY);
  }
  return *this;
}

template <> irtkGenericImage<irtkGreyPixel>& irtkGenericImage<irtkGreyPixel>::operator=
(const irtkGenericImage<irtkVector3D<double> >& image)
{
  int i, n;
  irtkGreyPixel* ptr1;
  irtkVector3D<double>* ptr2, min, max;
  double lmin, lmax;

  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(min._x*min._x + min._y*min._y + min._z*min._z);
  lmax = sqrt(max._x*max._x + max._y*max._y + max._z*max._z);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z);
    ptr1[i] = round((l - lmin)/(double)(lmax - lmin)*(MAX_GREY - MIN_GREY) + MIN_GREY);
  }
  return *this;
}

template <> irtkGenericImage<irtkRealPixel>& irtkGenericImage<irtkRealPixel>::operator=
(const irtkGenericImage<irtkBytePixel> &image)
{
  int i, n;
  irtkRealPixel *ptr1;
  irtkBytePixel *ptr2;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkRealPixel>& irtkGenericImage<irtkRealPixel>::operator=
(const irtkGenericImage<irtkGreyPixel> &image)
{
  int i, n;
  irtkRealPixel *ptr1;
  irtkGreyPixel *ptr2;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkRealPixel>& irtkGenericImage<irtkRealPixel>::operator=
(const irtkGenericImage<irtkRealPixel> &image)
{
  int i, n;
  irtkRealPixel *ptr1, *ptr2;

  if (this == &image) return *this;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkRealPixel>& irtkGenericImage<irtkRealPixel>::operator=
(const irtkGenericImage<irtkVector3D<char> >& image)
{
  int i, n;
  irtkRealPixel *ptr1;
  irtkVector3D<char> *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = static_cast<irtkRealPixel>(l);
  }
  return *this;
}

template <> irtkGenericImage<irtkRealPixel>& irtkGenericImage<irtkRealPixel>::operator=
(const irtkGenericImage<irtkVector3D<short> >& image)
{
  int i, n;
  irtkRealPixel *ptr1;
  irtkVector3D<short> *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = static_cast<irtkRealPixel>(l);
  }
  return *this;
}

template <> irtkGenericImage<irtkRealPixel>& irtkGenericImage<irtkRealPixel>::operator=
(const irtkGenericImage<irtkVector3D<float> >& image)
{
  int i, n;
  irtkRealPixel *ptr1;
  irtkVector3D<float> *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z);
    ptr1[i] = static_cast<irtkRealPixel>(l);
  }
  return *this;
}

template <> irtkGenericImage<irtkRealPixel>& irtkGenericImage<irtkRealPixel>::operator=
(const irtkGenericImage<irtkVector3D<double> >& image)
{
  int i, n;
  irtkRealPixel *ptr1;
  irtkVector3D<double> *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z);
    ptr1[i] = static_cast<irtkRealPixel>(l);
  }
  return *this;
}

template <> irtkGenericImage<irtkRGBPixel>& irtkGenericImage<irtkRGBPixel>::operator=
(const irtkGenericImage<irtkBytePixel> &image)
{
  int i, n;
  irtkRGBPixel *ptr1;
  irtkBytePixel *ptr2;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._r = ptr2[i];
    ptr1[i]._g = ptr2[i];
    ptr1[i]._b = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkRGBPixel>& irtkGenericImage<irtkRGBPixel>::operator=
(const irtkGenericImage<irtkGreyPixel> &image)
{
  int i, n;
  irtkRGBPixel *ptr1;
  irtkGreyPixel *ptr2;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if ((ptr2[i] >= 0) && (ptr2[i] < 256)) {
      ptr1[i]._r = ptr2[i];
      ptr1[i]._g = ptr2[i];
      ptr1[i]._b = ptr2[i];
    } else {
      if (ptr2[i] < 0) {
        ptr1[i]._r = 0;
        ptr1[i]._g = 0;
        ptr1[i]._b = 0;
      } else {
        ptr1[i]._r = 255;
        ptr1[i]._g = 255;
        ptr1[i]._b = 255;
      }
    }
  }
  return *this;
}

template <> irtkGenericImage<irtkRGBPixel>& irtkGenericImage<irtkRGBPixel>::operator=
(const irtkGenericImage<irtkRealPixel> &image)
{
  int i, n;
  irtkRGBPixel *ptr1;
  irtkRealPixel *ptr2;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if ((ptr2[i] >= 0) && (ptr2[i] < 256)) {
      ptr1[i]._r = round(ptr2[i]);
      ptr1[i]._g = round(ptr2[i]);
      ptr1[i]._b = round(ptr2[i]);
    } else {
      if (ptr2[i] < 0) {
        ptr1[i]._r = 0;
        ptr1[i]._g = 0;
        ptr1[i]._b = 0;
      } else {
        ptr1[i]._r = 255;
        ptr1[i]._g = 255;
        ptr1[i]._b = 255;
      }
    }
  }
  return *this;
}

template <> irtkGenericImage<irtkRGBPixel>& irtkGenericImage<irtkRGBPixel>::operator=
(const irtkGenericImage<irtkRGBPixel> &image)
{
  int i, n;
  irtkRGBPixel *ptr1, *ptr2;

  if (this == &image) return *this;

  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<char> >& irtkGenericImage<irtkVector3D<char> >::operator=
(const irtkGenericImage<irtkBytePixel>& image)
{
  int i, n;
  irtkVector3D<char> *ptr1;
  irtkBytePixel *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<short> >& irtkGenericImage<irtkVector3D<short> >::operator=
(const irtkGenericImage<irtkBytePixel>& image)
{
  int i, n;
  irtkVector3D<short> *ptr1;
  irtkBytePixel *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<float> >& irtkGenericImage<irtkVector3D<float> >::operator=
(const irtkGenericImage<irtkBytePixel>& image)
{
  int i, n;
  irtkVector3D<float> *ptr1;
  irtkBytePixel *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<double> >& irtkGenericImage<irtkVector3D<double> >::operator=
(const irtkGenericImage<irtkBytePixel>& image)
{
  int i, n;
  irtkVector3D<double> *ptr1;
  irtkBytePixel *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<char> >& irtkGenericImage<irtkVector3D<char> >::operator=
(const irtkGenericImage<irtkGreyPixel>& image)
{
  int i, n;
  irtkVector3D<char> *ptr1;
  irtkGreyPixel *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = static_cast<char>(ptr2[i]);
    ptr1[i]._y = static_cast<char>(ptr2[i]);
    ptr1[i]._z = static_cast<char>(ptr2[i]);
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<short> >& irtkGenericImage<irtkVector3D<short> >::operator=
(const irtkGenericImage<irtkGreyPixel>& image)
{
  int i, n;
  irtkVector3D<short> *ptr1;
  irtkGreyPixel *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<float> >& irtkGenericImage<irtkVector3D<float> >::operator=
(const irtkGenericImage<irtkGreyPixel>& image)
{
  int i, n;
  irtkVector3D<float> *ptr1;
  irtkGreyPixel *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<double> >& irtkGenericImage<irtkVector3D<double> >::operator=
(const irtkGenericImage<irtkGreyPixel>& image)
{
  int i, n;
  irtkVector3D<double> *ptr1;
  irtkGreyPixel *ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<char> >& irtkGenericImage<irtkVector3D<char> >::operator=
(const irtkGenericImage<irtkRealPixel>& image)
{
  int i, n;
  irtkVector3D<char>* ptr1;
  irtkRealPixel* ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = static_cast<char>(ptr2[i]);
    ptr1[i]._y = static_cast<char>(ptr2[i]);
    ptr1[i]._z = static_cast<char>(ptr2[i]);
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<short> >& irtkGenericImage<irtkVector3D<short> >::operator=
(const irtkGenericImage<irtkRealPixel>& image)
{
  int i, n;
  irtkVector3D<short>* ptr1;
  irtkRealPixel* ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = static_cast<short>(ptr2[i]);
    ptr1[i]._y = static_cast<short>(ptr2[i]);
    ptr1[i]._z = static_cast<short>(ptr2[i]);
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<float> >& irtkGenericImage<irtkVector3D<float> >::operator=
(const irtkGenericImage<irtkRealPixel>& image)
{
  int i, n;
  irtkVector3D<float>* ptr1;
  irtkRealPixel* ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<double> >& irtkGenericImage<irtkVector3D<double> >::operator=
(const irtkGenericImage<irtkRealPixel>& image)
{
  int i, n;
  irtkVector3D<double>* ptr1;
  irtkRealPixel* ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<char> >& irtkGenericImage<irtkVector3D<char> >::operator=
(const irtkGenericImage<irtkVector3D<char> >& image)
{
  int i, n;
  irtkVector3D<char>* ptr1;
  irtkVector3D<char>* ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<short> >& irtkGenericImage<irtkVector3D<short> >::operator=
(const irtkGenericImage<irtkVector3D<short> >& image)
{
  int i, n;
  irtkVector3D<short>* ptr1;
  irtkVector3D<short>* ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<float> >& irtkGenericImage<irtkVector3D<float> >::operator=
(const irtkGenericImage<irtkVector3D<float> >& image)
{
  int i, n;
  irtkVector3D<float>* ptr1;
  irtkVector3D<float>* ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkVector3D<double> >& irtkGenericImage<irtkVector3D<double> >::operator=
(const irtkGenericImage<irtkVector3D<double> >& image)
{
  int i, n;
  irtkVector3D<double>* ptr1;
  irtkVector3D<double>* ptr2;

  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <> irtkGenericImage<irtkBytePixel>::irtkGenericImage(const irtkGenericImage<irtkBytePixel> &image)
{
  int i, n;
  irtkBytePixel *ptr1, *ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <> irtkGenericImage<irtkBytePixel>::irtkGenericImage(const irtkGenericImage<irtkGreyPixel> &image)
{
  int i, n;
  irtkBytePixel *ptr1;
  irtkGreyPixel *ptr2, min, max;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = round(MAX_BYTE*(ptr2[i] - min)/(double)(max - min));
  }
}

template <> irtkGenericImage<irtkBytePixel>::irtkGenericImage(const irtkGenericImage<irtkRealPixel> &image)
{
  int i, n;
  irtkBytePixel *ptr1;
  irtkRealPixel *ptr2, min, max;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = round(MAX_BYTE*(ptr2[i] - min)/(double)(max - min));
  }
}

template <> irtkGenericImage<irtkBytePixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<char> >& image)
{
  int i, n;
  irtkBytePixel* ptr1;
  irtkVector3D<char>* ptr2, min, max;
  double lmin, lmax;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round(MAX_BYTE*(l - lmin)/(lmax - lmin));
  }
}

template <> irtkGenericImage<irtkBytePixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<short> >& image)
{
  int i, n;
  irtkBytePixel* ptr1;
  irtkVector3D<short>* ptr2, min, max;
  double lmin, lmax;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round(MAX_BYTE*(l - lmin)/(lmax - lmin));
  }
}

template <> irtkGenericImage<irtkBytePixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<float> >& image)
{
  int i, n;
  irtkBytePixel* ptr1;
  irtkVector3D<float>* ptr2, min, max;
  double lmin, lmax;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round(MAX_BYTE*(l - lmin)/(lmax - lmin));
  }
}

template <> irtkGenericImage<irtkBytePixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<double> >& image)
{
  int i, n;
  irtkBytePixel* ptr1;
  irtkVector3D<double>* ptr2, min, max;
  double lmin, lmax;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(min._x*min._x + min._y*min._y + min._z*min._z);
  lmax = sqrt(max._x*max._x + max._y*max._y + max._z*max._z);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z);
    ptr1[i] = round(MAX_BYTE*(l - lmin)/(lmax - lmin));
  }
}

template <> irtkGenericImage<irtkGreyPixel>::irtkGenericImage(const irtkGenericImage<irtkBytePixel> &image)
{
  int i, n;
  irtkGreyPixel *ptr1;
  irtkBytePixel *ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <> irtkGenericImage<irtkGreyPixel>::irtkGenericImage(const irtkGenericImage<irtkGreyPixel> &image)
{
  int i, n;
  irtkGreyPixel *ptr1, *ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <> irtkGenericImage<irtkGreyPixel>::irtkGenericImage(const irtkGenericImage<irtkRealPixel> &image)
{
  int i, n;
  irtkGreyPixel *ptr1;
  irtkRealPixel *ptr2, min, max;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = round(((ptr2[i] - min)/(double)(max - min)) *
                    (MAX_GREY - MIN_GREY) + MIN_GREY);
  }
}

template <> irtkGenericImage<irtkGreyPixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<char> >& image)
{
  int i, n;
  irtkGreyPixel* ptr1;
  irtkVector3D<char>* ptr2, min, max;
  double lmin, lmax;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round((l - lmin)/(lmax - lmin)*(MAX_GREY - MIN_GREY) + MIN_GREY);
  }
}

template <> irtkGenericImage<irtkGreyPixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<short> >& image)
{
  int i, n;
  irtkGreyPixel* ptr1;
  irtkVector3D<short>* ptr2, min, max;
  double lmin, lmax;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(static_cast<double>(min._x*min._x + min._y*min._y + min._z*min._z));
  lmax = sqrt(static_cast<double>(max._x*max._x + max._y*max._y + max._z*max._z));
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z));
    ptr1[i] = round((l - lmin)/(lmax - lmin)*(MAX_GREY - MIN_GREY) + MIN_GREY);
  }
}

template <> irtkGenericImage<irtkGreyPixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<float> >& image)
{
  int i, n;
  irtkGreyPixel* ptr1;
  irtkVector3D<float>* ptr2, min, max;
  double lmin, lmax;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(min._x*min._x + min._y*min._y + min._z*min._z);
  lmax = sqrt(max._x*max._x + max._y*max._y + max._z*max._z);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z);
    ptr1[i] = round((l - lmin)/(lmax - lmin)*(MAX_GREY - MIN_GREY) + MIN_GREY);
  }
}

template <> irtkGenericImage<irtkGreyPixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<double> >& image)
{
  int i, n;
  irtkGreyPixel* ptr1;
  irtkVector3D<double>* ptr2, min, max;
  double lmin, lmax;

  _matrix = NULL;
  this->Initialize(image);
  image.GetMinMax(&min, &max);
  lmin = sqrt(min._x*min._x + min._y*min._y + min._z*min._z);
  lmax = sqrt(max._x*max._x + max._y*max._y + max._z*max._z);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    double l = sqrt(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z);
    ptr1[i] = round((l - lmin)/(lmax - lmin)*(MAX_GREY - MIN_GREY) + MIN_GREY);
  }
}

template <> irtkGenericImage<irtkRealPixel>::irtkGenericImage(const irtkGenericImage<irtkBytePixel> &image)
{
  int i, n;
  irtkRealPixel *ptr1;
  irtkBytePixel *ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <> irtkGenericImage<irtkRealPixel>::irtkGenericImage(const irtkGenericImage<irtkGreyPixel> &image)
{
  int i, n;
  irtkRealPixel *ptr1;
  irtkGreyPixel *ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <> irtkGenericImage<irtkRealPixel>::irtkGenericImage(const irtkGenericImage<irtkRealPixel> &image)
{
  int i, n;
  irtkRealPixel *ptr1, *ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <> irtkGenericImage<irtkRealPixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<char> >& image)
{
  int i, n;
  irtkRealPixel* ptr1;
  irtkVector3D<char>* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = static_cast<irtkRealPixel>(sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z)));
  }
}

template <> irtkGenericImage<irtkRealPixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<short> >& image)
{
  int i, n;
  irtkRealPixel* ptr1;
  irtkVector3D<short>* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = static_cast<irtkRealPixel>(sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z)));
  }
}

template <> irtkGenericImage<irtkRealPixel>::irtkGenericImage(const irtkGenericImage<irtkVector3D<float> >& image)
{
  int i, n;
  irtkRealPixel* ptr1;
  irtkVector3D<float>* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = static_cast<irtkRealPixel>(sqrt(static_cast<double>(ptr2[i]._x*ptr2[i]._x + ptr2[i]._y*ptr2[i]._y + ptr2[i]._z*ptr2[i]._z)));
  }
}

template <> irtkGenericImage<irtkRGBPixel >::irtkGenericImage(const irtkGenericImage<irtkBytePixel>& image)
{
  int i, n;
  irtkRGBPixel*  ptr1;
  irtkBytePixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._r = ptr2[i];
    ptr1[i]._g = ptr2[i];
    ptr1[i]._b = ptr2[i];
  }
}

template <> irtkGenericImage<irtkRGBPixel >::irtkGenericImage(const irtkGenericImage<irtkGreyPixel>& image)
{
  int i, n;
  irtkRGBPixel*  ptr1;
  irtkGreyPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._r = ptr2[i];
    ptr1[i]._g = ptr2[i];
    ptr1[i]._b = ptr2[i];
  }
}

template <> irtkGenericImage<irtkRGBPixel>::irtkGenericImage(const irtkGenericImage<irtkRealPixel>& image)
{
  int i, n;
  irtkRGBPixel*  ptr1;
  irtkRealPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if ((round(ptr2[i]) >= 0) && (round(ptr2[i]) < 256)) {
      ptr1[i]._r = round(ptr2[i]);
      ptr1[i]._g = round(ptr2[i]);
      ptr1[i]._b = round(ptr2[i]);
    } else {
      if (ptr2[i] < 0) {
        ptr1[i]._r = 0;
        ptr1[i]._g = 0;
        ptr1[i]._b = 0;
      } else {
        ptr1[i]._r = 255;
        ptr1[i]._g = 255;
        ptr1[i]._b = 255;
      }
    }
  }
}

template <> irtkGenericImage<irtkRGBPixel >::irtkGenericImage(const irtkGenericImage<irtkRGBPixel>& image)
{
  int i, n;
  irtkRGBPixel *ptr1, *ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._r = ptr2[i]._r;
    ptr1[i]._g = ptr2[i]._g;
    ptr1[i]._b = ptr2[i]._b;
  }
}

template <> irtkGenericImage<irtkVector3D<char> >::irtkGenericImage(const irtkGenericImage<irtkBytePixel>& image)
{
  int i, n;
  irtkVector3D<char>* ptr1;
  irtkBytePixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<short> >::irtkGenericImage(const irtkGenericImage<irtkBytePixel>& image)
{
  int i, n;
  irtkVector3D<short>* ptr1;
  irtkBytePixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<float> >::irtkGenericImage(const irtkGenericImage<irtkBytePixel>& image)
{
  int i, n;
  irtkVector3D<float>* ptr1;
  irtkBytePixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<double> >::irtkGenericImage(const irtkGenericImage<irtkBytePixel>& image)
{
  int i, n;
  irtkVector3D<double>* ptr1;
  irtkBytePixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<char> >::irtkGenericImage(const irtkGenericImage<irtkGreyPixel>& image)
{
  int i, n;
  irtkVector3D<char>* ptr1;
  irtkGreyPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = static_cast<char>(ptr2[i]);
    ptr1[i]._y = static_cast<char>(ptr2[i]);
    ptr1[i]._z = static_cast<char>(ptr2[i]);
  }
}

template <> irtkGenericImage<irtkVector3D<short> >::irtkGenericImage(const irtkGenericImage<irtkGreyPixel>& image)
{
  int i, n;
  irtkVector3D<short>* ptr1;
  irtkGreyPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<float> >::irtkGenericImage(const irtkGenericImage<irtkGreyPixel>& image)
{
  int i, n;
  irtkVector3D<float>* ptr1;
  irtkGreyPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<double> >::irtkGenericImage(const irtkGenericImage<irtkGreyPixel>& image)
{
  int i, n;
  irtkVector3D<double>* ptr1;
  irtkGreyPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<char> >::irtkGenericImage(const irtkGenericImage<irtkRealPixel>& image)
{
  int i, n;
  irtkVector3D<char>* ptr1;
  irtkRealPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = static_cast<char>(ptr2[i]);
    ptr1[i]._y = static_cast<char>(ptr2[i]);
    ptr1[i]._z = static_cast<char>(ptr2[i]);
  }
}

template <> irtkGenericImage<irtkVector3D<short> >::irtkGenericImage(const irtkGenericImage<irtkRealPixel>& image)
{
  int i, n;
  irtkVector3D<short>* ptr1;
  irtkRealPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = static_cast<short>(ptr2[i]);
    ptr1[i]._y = static_cast<short>(ptr2[i]);
    ptr1[i]._z = static_cast<short>(ptr2[i]);
  }
}

template <> irtkGenericImage<irtkVector3D<float> >::irtkGenericImage(const irtkGenericImage<irtkRealPixel>& image)
{
  int i, n;
  irtkVector3D<float>* ptr1;
  irtkRealPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = static_cast<float>(ptr2[i]);
    ptr1[i]._y = static_cast<float>(ptr2[i]);
    ptr1[i]._z = static_cast<float>(ptr2[i]);
  }
}

template <> irtkGenericImage<irtkVector3D<double> >::irtkGenericImage(const irtkGenericImage<irtkRealPixel>& image)
{
  int i, n;
  irtkVector3D<double>* ptr1;
  irtkRealPixel* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i]._x = ptr2[i];
    ptr1[i]._y = ptr2[i];
    ptr1[i]._z = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<char> >::irtkGenericImage(const irtkGenericImage<irtkVector3D<char> >& image)
{
  int i, n;
  irtkVector3D<char>* ptr1;
  irtkVector3D<char>* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<short> >::irtkGenericImage(const irtkGenericImage<irtkVector3D<short> >& image)
{
  int i, n;
  irtkVector3D<short>* ptr1;
  irtkVector3D<short>* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<float> >::irtkGenericImage(const irtkGenericImage<irtkVector3D<float> >& image)
{
  int i, n;
  irtkVector3D<float>* ptr1;
  irtkVector3D<float>* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <> irtkGenericImage<irtkVector3D<double> >::irtkGenericImage(const irtkGenericImage<irtkVector3D<double> >& image)
{
  int i, n;
  irtkVector3D<double>* ptr1;
  irtkVector3D<double>* ptr2;

  _matrix = NULL;
  this->Initialize(image);
  n = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}



#ifdef HAS_VTK

template <> int irtkGenericImage<irtkBytePixel>::ImageToVTKScalarType()
{
  return VTK_UNSIGNED_CHAR;
}

template <> int irtkGenericImage<irtkGreyPixel>::ImageToVTKScalarType()
{
  return VTK_SHORT;
}

template <> int irtkGenericImage<irtkRealPixel>::ImageToVTKScalarType()
{
  return VTK_FLOAT;
}

template <> int irtkGenericImage<irtkRGBPixel>::ImageToVTKScalarType()
{
  cerr << "irtkGenericImage<irtkRGBPixel>::ImageToVTKScalarType: Not a scalar" << endl;
  exit(1);
}

template <> int irtkGenericImage<irtkVector3D<char> >::ImageToVTKScalarType()
{
  std::cerr << "irtkGenericImage<irtkVector3D<char> >::ImageToVTKScalarType: Not a scalar." << std::endl;
  exit(1);
}

template <> int irtkGenericImage<irtkVector3D<short> >::ImageToVTKScalarType()
{
  std::cerr << "irtkGenericImage<irtkVector3D<short> >::ImageToVTKScalarType: Not a scalar." << std::endl;
  exit(1);
}

template <> int irtkGenericImage<irtkVector3D<float> >::ImageToVTKScalarType()
{
  std::cerr << "irtkGenericImage<irtkVector3D<float> >::ImageToVTKScalarType: Not a scalar." << std::endl;
  exit(1);
}

template <> int irtkGenericImage<irtkVector3D<double> >::ImageToVTKScalarType()
{
  std::cerr << "irtkGenericImage<irtkVector3D<double> >::ImageToVTKScalarType: Not a scalar." << std::endl;
  exit(1);
}

template <class Type> void irtkGenericImage<Type>::ImageToVTK(vtkStructuredPoints *vtk)
{
  int i, n;
  double x, y, z;
  Type *ptr1, *ptr2;

  // Calculate the VTK origin of an IRTK image
  x = 0;
  y = 0;
  z = 0;
  this->ImageToWorld(x, y, z);

  // Allocate the VTK image
  vtk->SetOrigin(x, y, z);
  vtk->SetDimensions(_x, _y, _z);
  vtk->SetSpacing(_dx, _dy, _dz);
  vtk->SetScalarType(this->ImageToVTKScalarType());
  vtk->AllocateScalars();

  // Initialize the VTK image
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = (Type *)vtk->GetScalarPointer();
  for (i = 0; i < this->GetNumberOfVoxels(); i++) {
    *ptr2 = *ptr1;
    ptr1++;
    ptr2++;
  }
}

template <class Type> void irtkGenericImage<Type>::VTKToImage(vtkStructuredPoints *vtk)
{}

#endif

#include <irtkTemplate.h>
