/*=========================================================================
 
  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details
 
=========================================================================*/

#define _IMPLEMENTS_GENERICIMAGE_

#include <irtkImage.h>

#include <irtkFileToImage.h>
#include <irtkImageToFile.h>

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(void) : irtkBaseImage()
{
  _attr._x = 0;
  _attr._y = 0;
  _attr._z = 0;
  _attr._t = 0;

  // Initialize data
  _matrix  = NULL;
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, int t) : irtkBaseImage()
{
  irtkImageAttributes attr;

  attr._x = x;
  attr._y = y;
  attr._z = z;
  attr._t = t;

  // Initialize data
  _matrix = NULL;

  // Initialize rest of class
  this->Initialize(attr);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(char *filename)
{
  // Initialize data
  _matrix = NULL;

  // Read image
  this->Read(filename);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(const irtkImageAttributes &attr) : irtkBaseImage()
{
  // Initialize data
  _matrix  = NULL;

  // Initialize rest of class
  this->Initialize(attr);
}

template <class VoxelType> irtkGenericImage<VoxelType>::irtkGenericImage(const irtkGenericImage &image) : irtkBaseImage()
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  // Initialize data
  _matrix = NULL;

  // Initialize rest of class
  this->Initialize(image._attr);

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
}

template <class VoxelType> template <class VoxelType2> irtkGenericImage<VoxelType>::irtkGenericImage(const irtkGenericImage<VoxelType2> &image)
{
  int i, n;
  VoxelType  *ptr1;
  VoxelType2 *ptr2;

  // Initialize data
  _matrix = NULL;

  // Initialize rest of class
  this->Initialize(image.GetImageAttributes());

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = static_cast<VoxelType>(ptr2[i]);
  }
}

template <class VoxelType> irtkGenericImage<VoxelType>::~irtkGenericImage(void)
{
  if (_matrix != NULL) {
    Deallocate<VoxelType>(_matrix);
    _matrix = NULL;
  }
  _attr._x = 0;
  _attr._y = 0;
  _attr._z = 0;
  _attr._t = 0;
}

template <> const char *irtkGenericImage<char>::NameOfClass()
{
  return "irtkGenericImage<char>";
}

template <> const char *irtkGenericImage<unsigned char>::NameOfClass()
{
  return "irtkGenericImage<unsigned char>";
}

template <> const char *irtkGenericImage<short>::NameOfClass()
{
  return "irtkGenericImage<short>";
}

template <> const char *irtkGenericImage<unsigned short>::NameOfClass()
{
  return "irtkGenericImage<unsigned short>";
}

template <> const char *irtkGenericImage<int>::NameOfClass()
{
  return "irtkGenericImage<int>";
}

template <> const char *irtkGenericImage<unsigned int>::NameOfClass()
{
  return "irtkGenericImage<unsigned int>";
}

template <> const char *irtkGenericImage<float>::NameOfClass()
{
  return "irtkGenericImage<float>";
}

template <> const char *irtkGenericImage<double>::NameOfClass()
{
  return "irtkGenericImage<double>";
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize(const irtkImageAttributes &attr)
{
  // Free memory
  if ((_attr._x != attr._x) || (_attr._y != attr._y) || (_attr._z != attr._z) || (_attr._t != attr._t)) {
    // Free old memory
    if (_matrix != NULL) Deallocate<VoxelType>(_matrix);
    // Allocate new memory
    if (attr._x*attr._y*attr._z*attr._t > 0) {
      _matrix = Allocate(_matrix, attr._x, attr._y, attr._z, attr._t);
    } else {
      _matrix = NULL;
    }
  }

  // Initialize base class
  this->irtkBaseImage::Update(attr);

  // Initialize voxels
  *this = VoxelType();
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Clear()
{
	// Free memory
	if (_matrix != NULL)
		Deallocate<VoxelType>(_matrix);

  _attr._x = 0;
  _attr._y = 0;
  _attr._z = 0;
  _attr._t = 0;

}

template <class VoxelType> void irtkGenericImage<VoxelType>::Read(const char *filename)
{
  irtkBaseImage *image;

  // Allocate file reader
  irtkFileToImage *reader = irtkFileToImage::New(filename);

  // Get output
  image = reader->GetOutput();

  // Convert image
  switch (reader->GetDataType()) {
          
  case IRTK_VOXEL_CHAR: { *this = *(dynamic_cast<irtkGenericImage<char> *>(image)); } break;

  case IRTK_VOXEL_UNSIGNED_CHAR: { *this = *(dynamic_cast<irtkGenericImage<unsigned char> *>(image)); } break;

  case IRTK_VOXEL_SHORT: { *this = *(dynamic_cast<irtkGenericImage<short> *>(image)); } break;

  case IRTK_VOXEL_UNSIGNED_SHORT: { *this = *(dynamic_cast<irtkGenericImage<unsigned short> *>(image)); } break;
    
  case IRTK_VOXEL_INT: { *this = *(dynamic_cast<irtkGenericImage<int> *>(image)); } break;

  case IRTK_VOXEL_UNSIGNED_INT: { *this = *(dynamic_cast<irtkGenericImage<unsigned int> *>(image)); } break;

  case IRTK_VOXEL_FLOAT: { *this = *(dynamic_cast<irtkGenericImage<float> *>(image)); } break;

  case IRTK_VOXEL_DOUBLE: { *this = *(dynamic_cast<irtkGenericImage<double> *>(image)); } break;

  default:
      cout << "irtkGenericImage::GetOutput: Unknown voxel type" << endl;
  }

  if (reader->GetSlope() != 0) {
      switch (this->GetScalarType()) {

      case IRTK_VOXEL_FLOAT: {
          *this *= static_cast<float>(reader->GetSlope());
          *this += static_cast<float>(reader->GetIntercept());
      }
          break;

      case IRTK_VOXEL_DOUBLE: {
          *this *= static_cast<double>(reader->GetSlope());
          *this += static_cast<double>(reader->GetIntercept());
      }
          break;

      default:
          if ((reader->GetSlope() != 1) || (reader->GetIntercept() != 0)) {
              cerr << this->NameOfClass() << "::Read: Ignore slope and intercept, use irtkGenericImage<float> or " << endl;
              cerr << "irtkGenericImage<double> instead" << endl;
          }
      }
  }
  
  // Delete reader
  delete reader;

  // Delete image
  delete image;
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Write(const char *filename)
{
  // Allocate file reader
  irtkImageToFile *writer = irtkImageToFile::New(filename);

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
      if (ptr[i] < *min)
        *min = ptr[i];
      if (ptr[i] > *max)
        *max = ptr[i];
    }
  }
}

template <class VoxelType> VoxelType irtkGenericImage<VoxelType>::GetAverage(int toggle) const
{
  float average = 0;
  int i, n, m;
  VoxelType *ptr;

  // Initialize pixels
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();

  m = 0;
  if (n > 0) {
    for (i = 0; i < n; i++) {
      if(*ptr > 0) {
        m ++;
      }
      ptr ++;
    }
    ptr = this->GetPointerToVoxels();
    for (i = 0; i < n; i++) {
      if(toggle == 1){
		  if(*ptr > 0) {
			  average += (float)((VoxelType)*ptr)/(float)m;
		  }
      }
      else{
        average += (float)((VoxelType)*ptr)/(float)n;
      }
      ptr ++;
    }
  }
  return average;
}

template <class VoxelType> VoxelType irtkGenericImage<VoxelType>::GetSD(int toggle) const
{
  // Initialize pixels
  float average = 0, std = 0;
  int i, n;
  VoxelType *ptr;
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  average = this->GetAverage(toggle);

  if (n > 0) {
    for (i = 0; i < n; i++) {
      if(*ptr > 0 && toggle == 1) {
        std += pow((double)*ptr - average,2)/(double)n;
      } else {
        std += pow((double)*ptr - average,2)/(double)n;
      }
      ptr ++;
    }
  }
  return sqrt(std);
}

template <class VoxelType> VoxelType irtkGenericImage<VoxelType>::GetAveragePad(double padding) const
{
  VoxelType average = 0;
  int i, n, m;
  VoxelType *ptr;

  // Initialize pixels
  n = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();

  m = 0;
  for (i = 0; i < n; i++) {
	  if(*ptr > padding) {
		  average += (float)((VoxelType)*ptr);
		  m++;
	  }
      ptr++;
  }

  if (m > 0) average /= double(m);

  return average;
}

template <class VoxelType> VoxelType irtkGenericImage<VoxelType>::GetSDPad(double padding) const
{
  // Initialize pixels
  VoxelType average = 0, std = 0;
  int i, n, m;
  VoxelType *ptr;
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  average = this->GetAveragePad(padding);

  m = 0;
  for (i = 0; i < n; i++) {
      if(*ptr > padding) {
        std += pow((double)*ptr - average,2);
        m++;
      }
      ptr++;
  }

  if (m > 0) std /= double(m);

  return sqrt((double)std);
}

template <class VoxelType> void irtkGenericImage<VoxelType>::GetMaxPosition(irtkPoint& p, int ds, int) const
{
  double i, j, k;
  VoxelType *ptr;
  double x, y , z;

  this->WorldToImage(p);
  VoxelType max = 0;
  x = round(p._x); y = round(p._y); z = round(p._z);
  k = round(p._z);
  ptr = this->GetPointerToVoxels();
  for (j = round(p._y) - ds; j < round(p._y) +ds + 1; j++) {
    for (i = round(p._x) - ds; i < round(p._x)+ds + 1; i++) {
      // Initialize pixels
      if (max < *ptr) {
        max = *ptr;
        x = i; y = j; z = k;
      }
      ptr++;
    }
  }
  p._x = x; p._y = y; p._z = z;
  this->ImageToWorld(p);
}

template <class VoxelType> void irtkGenericImage<VoxelType>::GravityCenter(irtkPoint& p, int ds, int t) const
{
  double i, j, k;
  //VoxelType *ptr = new VoxelType;
  VoxelType *ptr;
  double x,y,z;
  double si,sj,sk;
  double sweight;

  //irtkInterpolateImageFunction *interpolator = irtkInterpolateImageFunction::New(Interpolation_CSpline, (irtkBaseImage*)this);
  // Setup interpolation for the source image
  //interpolator->SetInput((irtkImage*)this);
  //interpolator->Initialize();

  this->WorldToImage(p);
  si = 0; sj = 0; sk = 0; sweight = 0;
  x = round(p._x); y = round(p._y); z = round(p._z);
  for (k = round(p._z) - ds; k < round(p._z) +ds + 1; k++) {
    for (j = round(p._y) - ds; j < round(p._y) +ds + 1; j++) {
      for (i = round(p._x) - ds; i < round(p._x)+ds + 1; i++) {
        // Initialize pixels
        ptr = this->GetPointerToVoxels(i,j,k,t);
        si += *ptr * i;
        sj += *ptr * j;
        sk += *ptr * k;
        sweight += *ptr;
      }
    }
  }
  x = si/sweight; y = sj/sweight; z = sk/sweight;
  p._x = x; p._y = y; p._z = z;
  this->ImageToWorld(p);
  //delete interpolator;
  //delete ptr;
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
      if (ptr[i]>pad) {
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

template <class VoxelType> void irtkGenericImage<VoxelType>::PutMinMax(VoxelType min, VoxelType max)
{
  int i, n;
  VoxelType *ptr, min_val, max_val;

  // Get lower and upper bound
  this->GetMinMax(&min_val, &max_val);

  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = VoxelType(((ptr[i] - min_val) / double(max_val - min_val)) *
                       (max - min) + min);
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Saturate( double q0, double q1 )
{
  if ((q0==0) && (q1==1))
    return;
  if ((q0<0) || (q1>1)) {
    stringstream msg;
    msg << "irtkGenericImage<VoxelType>::Saturate: Parameter out of range 0.0 - 1.0\n";
    cerr << msg.str();
    throw irtkException( msg.str(),
                         __FILE__,
                         __LINE__ );
  }
  
  int i, n;
  VoxelType *ptr, q0_val, q1_val;
  
  // find quantiles
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();  
  vector<VoxelType> voxel_tmp(n);
  for (i=0;i<n;i++)
    voxel_tmp[i] = ptr[i];

  sort(voxel_tmp.begin(),voxel_tmp.end());

  q0_val = voxel_tmp[round((voxel_tmp.size()-1)*q0)];
  q1_val = voxel_tmp[round((voxel_tmp.size()-1)*q1)];

  for (i = 0; i < n; i++) {
    if (ptr[i] < q0_val)
      ptr[i] = q0_val;
    if (ptr[i] > q1_val)
      ptr[i] = q1_val;
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::Quantiles( vector<double> &q, 
                                                                        vector<VoxelType> &res )
{
  int i, n;
  VoxelType *ptr;
  
  // find quantiles
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();  
  vector<VoxelType> voxel_tmp(n);
  for (i=0;i<n;i++)
    voxel_tmp[i] = ptr[i];

  sort(voxel_tmp.begin(),voxel_tmp.end());

  res.resize( q.size() );
  for (i=0;i<q.size();i++)
    res[i] = voxel_tmp[round((voxel_tmp.size()-1)*q[i])];
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::GetRegion(int k, int m) const
{
  int i, j;
  double x1, y1, z1, t1, x2, y2, z2, t2;

  if ((k < 0) || (k >= _attr._z) || (m < 0) || (m >= _attr._t)) {
      stringstream msg;
      msg << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
  }

  // Initialize
  irtkImageAttributes attr = this->_attr;
  attr._z = 1;
  attr._t = 1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  irtkGenericImage<VoxelType> image(attr);

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
  for (j = 0; j < _attr._y; j++) {
    for (i = 0; i < _attr._x; i++) {
      image._matrix[0][0][j][i] = _matrix[m][k][j][i];
    }
  }
  return image;
}

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::GetFrame(int l) const
{
  int i, j, k;

  if ((l < 0) || (l >= _attr._t)) {
      stringstream msg;
      msg << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
  }

  // Initialize
  irtkImageAttributes attr = this->_attr;
  attr._t = 1;
  irtkGenericImage<VoxelType> image(attr);

  // Copy region
  for (k = 0; k < _attr._z; k++) {
    for (j = 0; j < _attr._y; j++) {
      for (i = 0; i < _attr._x; i++) {
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
      (i2 > _attr._x) || (j2 > _attr._y) || (k2 > _attr._z)) {
      stringstream msg;
      msg << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
  }

  // Initialize
  irtkImageAttributes attr = this->_attr;
  attr._x = i2 - i1;
  attr._y = j2 - j1;
  attr._z = k2 - k1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  irtkGenericImage<VoxelType> image(attr);

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
  for (l = 0; l < _attr._t; l++) {
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
      (i2 > _attr._x) || (j2 > _attr._y) || (k2 > _attr._z) || (l2 > _attr._t)) {
      stringstream msg;
      msg << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
  }

  // Initialize
  irtkImageAttributes attr = this->_attr;
  attr._x = i2 - i1;
  attr._y = j2 - j1;
  attr._z = k2 - k1;
  attr._t = l2 - l1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  irtkGenericImage<VoxelType> image(attr);

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

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator=(const irtkGenericImage<VoxelType> &image)
{
  int i, n;
  VoxelType  *ptr1, *ptr2;

  if (this == &image) return *this;

  this->Initialize(image._attr);
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = ptr2[i];
  }
  return *this;
}

template <class VoxelType> template <class VoxelType2> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator=(const irtkGenericImage<VoxelType2> &image)
{
  int i, n;
  VoxelType  *ptr1;
  VoxelType2 *ptr2;

  this->Initialize(image.GetImageAttributes());
  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr1[i] = static_cast<VoxelType>(ptr2[i]);
  }
  return *this;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator+=(const irtkGenericImage<VoxelType> &image)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  if (!(this->GetImageAttributes() == image.GetImageAttributes())) {
      stringstream msg;
      msg << "irtkGenericImage<VoxelType>::operator+=: Size mismatch in images\n";
      cerr << msg.str();
      this->GetImageAttributes().Print();
      image.GetImageAttributes().Print();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
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

  if (!(this->GetImageAttributes() == image.GetImageAttributes())) {
      stringstream msg;
      msg << "irtkGenericImage<VoxelType>::operator-=: Size mismatch in images\n";
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
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

  if (!(this->GetImageAttributes() == image.GetImageAttributes())) {
      stringstream msg;
      msg << "irtkGenericImage<VoxelType>::operator*=: Size mismatch in images\n";
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
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

  if (!(this->GetImageAttributes() == image.GetImageAttributes())) {
      stringstream msg;
      msg << "irtkGenericImage<VoxelType>::operator/=: Size mismatch in images\n";
      cerr << msg.str();
      throw irtkException( msg.str(),
                           __FILE__,
                           __LINE__ );
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

template <class VoxelType> bool irtkGenericImage<VoxelType>::operator==(const irtkGenericImage<VoxelType> &image)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  if (!(this->GetImageAttributes() == image.GetImageAttributes())) {
    return false;
  }

  n    = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if (ptr1[i] != ptr2[i]) return false;
  }
  return true;
}

template <class VoxelType> irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator=(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr;
  
  n   = this->GetNumberOfVoxels();
  ptr = this->GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i] = pixel;
  }
  return *this;
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

template <class VoxelType> irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator!=(VoxelType pixel)
{
  int i, n;
  VoxelType *ptr1, *ptr2;

  irtkGenericImage<VoxelType> image(*this);

  n   = this->GetNumberOfVoxels();
  ptr1 = this->GetPointerToVoxels();
  ptr2 = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    if (ptr1[i] != pixel) {
      ptr2[i] = 1;
    }
    else {
        ptr2[i] = 0;
    }
  }
  return image;
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

  for (t = 0; t < _attr._t; t++) {
    for (z = 0; z < _attr._z; z++) {
      for (y = 0; y < _attr._y; y++) {
        for (x = 0; x < _attr._x / 2; x++) {
          swap(_matrix[t][z][y][x], _matrix[t][z][y][_attr._x-(x+1)]);
        }
      }
    }
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::ReflectY()
{
  int x, y, z, t;

  for (t = 0; t < _attr._t; t++) {
    for (z = 0; z < _attr._z; z++) {
      for (y = 0; y < _attr._y / 2; y++) {
        for (x = 0; x < _attr._x; x++) {
          swap(_matrix[t][z][y][x], _matrix[t][z][_attr._y-(y+1)][x]);
        }
      }
    }
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::ReflectZ()
{
  int x, y, z, t;

  for (t = 0; t < _attr._t; t++) {
    for (z = 0; z < _attr._z / 2; z++) {
      for (y = 0; y < _attr._y; y++) {
        for (x = 0; x < _attr._x; x++) {
          swap(_matrix[t][z][y][x], _matrix[t][_attr._z-(z+1)][y][x]);
        }
      }
    }
  }
}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipXY(int modifyOrigin)
{
  int i, j, k, m;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _attr._y, _attr._x, _attr._z, _attr._t);

  for (m = 0; m < _attr._t; m++) {
    for (k = 0; k < _attr._z; k++) {
      for (j = 0; j < _attr._y; j++) {
        for (i = 0; i < _attr._x; i++) {
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
  swap(_attr._x, _attr._y);

  // Swap voxel dimensions
  swap(_attr._dx, _attr._dy);

  if (modifyOrigin > 0){
  	// Swap origin coordinates
  	swap(_attr._xorigin, _attr._yorigin);
  }

  // Update transformation matrix
  this->UpdateMatrix();

}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipXZ(int modifyOrigin)
{
  int i, j, k, l;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _attr._z, _attr._y, _attr._x, _attr._t);

  for (l = 0; l < _attr._t; l++) {
    for (k = 0; k < _attr._z; k++) {
      for (j = 0; j < _attr._y; j++) {
        for (i = 0; i < _attr._x; i++) {
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
  swap(_attr._x, _attr._z);

  // Swap voxel dimensions
  swap(_attr._dx, _attr._dz);

  if (modifyOrigin > 0){
  	// Swap origin coordinates
  	swap(_attr._xorigin, _attr._zorigin);
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipYZ(int modifyOrigin)
{
  int i, j, k, l;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _attr._x, _attr._z, _attr._y, _attr._t);

  for (l = 0; l < _attr._t; l++) {
    for (k = 0; k < _attr._z; k++) {
      for (j = 0; j < _attr._y; j++) {
        for (i = 0; i < _attr._x; i++) {
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
  swap(_attr._y, _attr._z);

  // Swap voxel dimensions
  swap(_attr._dy, _attr._dz);

  if (modifyOrigin > 0){
  	// Swap origin coordinates
  	swap(_attr._yorigin, _attr._zorigin);
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipXT(int modifyOrigin)
{
  int i, j, k, m;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _attr._t, _attr._y, _attr._z, _attr._x);

  for (m = 0; m < _attr._t; m++) {
    for (k = 0; k < _attr._z; k++) {
      for (j = 0; j < _attr._y; j++) {
        for (i = 0; i < _attr._x; i++) {
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
  swap(_attr._x, _attr._t);

  // Swap voxel dimensions
  swap(_attr._dx, _attr._dt);

  if (modifyOrigin > 0){
  	// Swap origin coordinates
  	swap(_attr._xorigin, _attr._torigin);
  }

  // Update transformation matrix
  this->UpdateMatrix();

}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipYT(int modifyOrigin)
{
  int i, j, k, m;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _attr._x, _attr._t, _attr._z, _attr._y);

  for (m = 0; m < _attr._t; m++) {
    for (k = 0; k < _attr._z; k++) {
      for (j = 0; j < _attr._y; j++) {
        for (i = 0; i < _attr._x; i++) {
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
  swap(_attr._y, _attr._t);

  // Swap voxel dimensions
  swap(_attr._dy, _attr._dt);

  if (modifyOrigin > 0){
  	// Swap origin coordinates
  	swap(_attr._yorigin, _attr._torigin);
  }

  // Update transformation matrix
  this->UpdateMatrix();

}

template <class VoxelType> void irtkGenericImage<VoxelType>::FlipZT(int modifyOrigin)
{
  int i, j, k, m;
  VoxelType ****matrix = NULL;

  // Allocate memory
  matrix = Allocate(matrix, _attr._x, _attr._y, _attr._t, _attr._z);

  for (m = 0; m < _attr._t; m++) {
    for (k = 0; k < _attr._z; k++) {
      for (j = 0; j < _attr._y; j++) {
        for (i = 0; i < _attr._x; i++) {
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
  swap(_attr._z, _attr._t);

  // Swap voxel dimensions
  swap(_attr._dz, _attr._dt);

  if (modifyOrigin > 0){
  	// Swap origin coordinates
  	swap(_attr._zorigin, _attr._torigin);
  }

  // Update transformation matrix
  this->UpdateMatrix();

}

#ifdef HAS_VTK

template <> int irtkGenericImage<char>::ImageToVTKScalarType()
{
  return VTK_CHAR;
}

template <> int irtkGenericImage<unsigned char>::ImageToVTKScalarType()
{
  return VTK_UNSIGNED_CHAR;
}

template <> int irtkGenericImage<short>::ImageToVTKScalarType()
{
  return VTK_SHORT;
}

template <> int irtkGenericImage<unsigned short>::ImageToVTKScalarType()
{
  return VTK_UNSIGNED_SHORT;
}

template <> int irtkGenericImage<int>::ImageToVTKScalarType()
{
  return VTK_INT;
}

template <> int irtkGenericImage<unsigned int>::ImageToVTKScalarType()
{
  return VTK_UNSIGNED_INT;
}

template <> int irtkGenericImage<float>::ImageToVTKScalarType()
{
  return VTK_FLOAT;
}

template <> int irtkGenericImage<double>::ImageToVTKScalarType()
{
  return VTK_DOUBLE;
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
  vtk->SetDimensions(_attr._x, _attr._y, _attr._z);
  vtk->SetSpacing(_attr._dx, _attr._dy, _attr._dz);
#if VTK_MAJOR_VERSION >= 6
  vtkInformation* outInfo;
  vtk->SetScalarType(this->ImageToVTKScalarType(),outInfo);
  vtk->AllocateScalars(outInfo);
#else //VTK_MAJOR_VERSION >= 6
  vtk->SetScalarType(this->ImageToVTKScalarType());
  vtk->AllocateScalars();
#endif //VTK_MAJOR_VERSION >= 6

 

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

template <class Type> void irtkGenericImage<Type>::VTKToImage(vtkStructuredPoints *)
{}

#endif

#include <irtkTemplate.h>
