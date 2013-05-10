/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGENERICIMAGE_H

#define _IRTKGENERICIMAGE_H

#ifdef HAS_VTK

#include <vtkStructuredPoints.h>

#endif

/**
 * Generic class for 2D or 3D images
 *
 * This class implements generic 2D and 3D images. It provides functions
 * for accessing, reading, writing and manipulating images. This class can
 * be used for images with arbitrary voxel types using templates.
 */

template <class VoxelType> class irtkGenericImage : public irtkBaseImage
{

protected:

  /// Pointer to image data
  VoxelType ****_matrix;

public:

  /// Default constructor
  irtkGenericImage();

  /// Constructor from image file
  irtkGenericImage(char *);

  /// Constructor for given image size
  irtkGenericImage(int, int, int, int = 1);

  /// Copy constructor for image 
  irtkGenericImage(const irtkGenericImage &);

  /// Constructor for given image attributes
  irtkGenericImage(const irtkImageAttributes &);

  /// Copy constructor for image of different type
  template <class T> irtkGenericImage(const irtkGenericImage<T> &);

  /// Destructor
  ~irtkGenericImage(void);

  /// Initialize an image
  void Initialize(const irtkImageAttributes &);

  /// Clear an image
  void Clear();

  /// Read image from file
  void Read (const char *);

  /// Write image to file
  void Write(const char *);

  /// Minimum and maximum pixel values get accessor
  void GetMinMax(VoxelType *, VoxelType *) const;

  /// Average pixel values get accessor
  VoxelType GetAverage(int = 1) const;

  /// Standard Deviation of the pixels
  VoxelType GetSD(int = 1) const;

  /// Get Max Intensity position around the point
  void GetMaxPosition(irtkPoint &, int = 1, int = 0) const;

  /// Get Gravity center position of a given window
  void GravityCenter(irtkPoint &, int = 1, int = 0) const;

  /// Minimum and maximum pixel values get accessor with padding
  void GetMinMaxPad(VoxelType *, VoxelType *, VoxelType) const;

  /// Minimum and maximum pixel values put accessor
  void PutMinMax(VoxelType, VoxelType);

  /// Function for pixel access via pointers
  VoxelType *GetPointerToVoxels(int = 0, int = 0, int = 0, int = 0) const;

  /// Funnction to convert pixel to index
  int VoxelToIndex(int, int, int, int = 0) const;

  /// Function for pixel get access
  VoxelType   Get(int, int, int, int = 0) const;

  /// Function for pixel put access
  void   Put(int, int, int, VoxelType);

  /// Function for pixel put access
  void   Put(int, int, int, int, VoxelType);

  /// Function for pixel access from via operators
  VoxelType& operator()(int, int, int, int = 0);

  /// Function for image slice get access
  irtkGenericImage GetRegion(int z, int t) const;

  /// Function for image slice get access in certain region
  irtkGenericImage GetRegion(int x1, int y1, int z1, int x2, int y2, int z2) const;

  /// Function for image slice get access in certain region
  irtkGenericImage GetRegion(int x1, int y1, int z1, int t1, int x2, int y2, int z2, int t2) const;

  /// Function for image frame get access
  irtkGenericImage GetFrame(int t) const;

  //
  // Operators for image arithmetics
  //

  /// Copy operator for image
  irtkGenericImage<VoxelType>& operator= (const irtkGenericImage &);
  
  /// Copy operator for image
  template <class T> irtkGenericImage<VoxelType>& operator= (const irtkGenericImage<T> &);
  
  /// Addition operator
  irtkGenericImage  operator+ (const irtkGenericImage &);

  /// Addition operator (stores result)
  irtkGenericImage& operator+=(const irtkGenericImage &);

  /// Subtraction operator
  irtkGenericImage  operator- (const irtkGenericImage &);

  /// Subtraction operator (stores result)
  irtkGenericImage& operator-=(const irtkGenericImage &);

  /// Multiplication operator
  irtkGenericImage  operator* (const irtkGenericImage &);

  /// Multiplication operator (stores result)
  irtkGenericImage& operator*=(const irtkGenericImage &);

  /// Division operator
  irtkGenericImage  operator/ (const irtkGenericImage &);

  /// Division operator (stores result)
  irtkGenericImage& operator/=(const irtkGenericImage &);

  //
  // Operators for image and Type arithmetics
  //

  /// Set all pixels to a constant value
  irtkGenericImage& operator= (VoxelType);
  /// Addition operator for type
  irtkGenericImage  operator+ (VoxelType);
  /// Addition operator for type (stores result)
  irtkGenericImage& operator+=(VoxelType);
  /// Subtraction operator for type
  irtkGenericImage  operator- (VoxelType);
  /// Subtraction operator for type (stores result)
  irtkGenericImage& operator-=(VoxelType);
  /// Multiplication operator for type
  irtkGenericImage  operator* (VoxelType);
  /// Multiplication operator for type (stores result)
  irtkGenericImage& operator*=(VoxelType);
  /// Division operator for type
  irtkGenericImage  operator/ (VoxelType);
  /// Division operator for type (stores result)
  irtkGenericImage& operator/=(VoxelType);

  //
  // Operators for image thresholding
  //

  /// Threshold operator >  (sets all values >  given value to that value)
  irtkGenericImage  operator> (VoxelType);
  /// Threshold operator >= (sets all values >= given value to that value)
  irtkGenericImage& operator>=(VoxelType);
  /// Threshold operator <  (sets all values <  given value to that value)
  irtkGenericImage  operator< (VoxelType);
  /// Threshold operator <= (sets all values <= given value to that value)
  irtkGenericImage& operator<=(VoxelType);

  /// Comparison operators == (explicit negation yields != operator)
  bool operator==(const irtkGenericImage &);

  /// Comparison operator != (if _HAS_STL is defined, negate == operator)
  ///  bool operator!=(const irtkGenericImage &);

  irtkGenericImage operator!=(VoxelType);

  //
  // Reflections and axis flipping
  //

  /// Reflect image around x
  void ReflectX();
  /// Reflect image around y
  void ReflectY();
  /// Reflect image around z
  void ReflectZ();
  /// Flip x and y axis
  void FlipXY(int);
  /// Flip x and z axis
  void FlipXZ(int);
  /// Flip y and z axis
  void FlipYZ(int);
  /// Flip x and t axis
  void FlipXT(int);
  /// Flip y and t axis
  void FlipYT(int);
  /// Flip z and t axis
  void FlipZT(int);

  //
  // Conversions from and to VTK
  //

#ifdef HAS_VTK

  /// Return the VTK scalar image type of an IRTK image
  int ImageToVTKScalarType();

  /// Conversion to VTK structured points
  void ImageToVTK(vtkStructuredPoints *);

  /// Conversion from VTK structured points
  void VTKToImage(vtkStructuredPoints *);
#endif

  /// Function for pixel get access as double
  double GetAsDouble(int, int, int, int = 0) const;

  /// Function for pixel put access
  void   PutAsDouble(int, int, int, double);

  /// Function for pixel put access
  void   PutAsDouble(int, int, int, int, double);

  /// Returns the name of the image class
  const char *NameOfClass();
  
  /// Function for pixel access via pointers
  void *GetScalarPointer(int = 0, int = 0, int = 0, int = 0) const;

  /// Function which returns pixel scalar type
  virtual int GetScalarType() const;
  
  /// Function which returns the minimum value the pixel can hold without overflowing
  virtual double GetScalarTypeMin() const;

  /// Function which returns the minimum value the pixel can hold without overflowing
  virtual double GetScalarTypeMax() const;

};

template <class VoxelType> inline void irtkGenericImage<VoxelType>::Put(int x, int y, int z, VoxelType val)
{
#ifdef NO_BOUNDS
  _matrix[0][z][y][x] = val;
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0)) {
    cout << "irtkGenericImage<VoxelType>::Put: parameter out of range\n";
  } else {
    _matrix[0][z][y][x] = static_cast<VoxelType>(val);
  }
#endif
}

template <class VoxelType> inline void irtkGenericImage<VoxelType>::Put(int x, int y, int z, int t, VoxelType val)
{
#ifdef NO_BOUNDS
  _matrix[t][z][y][x] = val;
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0) || (t >= _attr._t) || (t < 0)) {
    cout << "irtkGenericImage<VoxelType>::Put: parameter out of range\n";
  } else {
    _matrix[t][z][y][x] = val;
  }
#endif
}

template <class VoxelType> inline void irtkGenericImage<VoxelType>::PutAsDouble(int x, int y, int z, double val)
{
  if (val > voxel_limits<VoxelType>::max()) val = voxel_limits<VoxelType>::max();
  if (val < voxel_limits<VoxelType>::min()) val = voxel_limits<VoxelType>::min();  

#ifdef NO_BOUNDS
  _matrix[0][z][y][x] = static_cast<VoxelType>(val);
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0) || (_attr._t > 0)) {
    cout << "irtkGenericImage<Type>::PutAsDouble: parameter out of range\n";
  } else {
    _matrix[0][z][y][x] = static_cast<VoxelType>(val);
  }
#endif
}

template <class VoxelType> inline void irtkGenericImage<VoxelType>::PutAsDouble(int x, int y, int z, int t, double val)
{
  if (val > voxel_limits<VoxelType>::max()) val = voxel_limits<VoxelType>::max();
  if (val < voxel_limits<VoxelType>::min()) val = voxel_limits<VoxelType>::min();  

#ifdef NO_BOUNDS
  _matrix[t][z][y][x] = static_cast<VoxelType>(val);
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0) || (t >= _attr._t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::PutAsDouble: parameter out of range\n";
  } else {
    _matrix[t][z][y][x] = static_cast<VoxelType>(val);
  }
#endif
}

template <class VoxelType> inline VoxelType irtkGenericImage<VoxelType>::Get(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return (_matrix[t][z][y][x]);
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0) || (t >= _attr._t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::Get: parameter out of range\n";
    return 0;
  } else {
    return(_matrix[t][z][y][x]);
  }
#endif
}

template <class VoxelType> inline double irtkGenericImage<VoxelType>::GetAsDouble(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return (static_cast<double>(_matrix[t][z][y][x]));
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0) || (t >= _attr._t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::GetAsDouble: parameter out of range\n";
    return 0;
  } else {
    return (static_cast<double>(_matrix[t][z][y][x]));
  }
#endif

}

template <class VoxelType> inline VoxelType& irtkGenericImage<VoxelType>::operator()(int x, int y, int z, int t)
{
#ifdef NO_BOUNDS
  return (_matrix[t][z][y][x]);
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0) || (t >= _attr._t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::(): parameter out of range\n";
    return _matrix[0][0][0][0];
  } else {
    return (_matrix[t][z][y][x]);
  }
#endif
}

template <class VoxelType> inline int irtkGenericImage<VoxelType>::VoxelToIndex(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return (&(_matrix[t][z][y][x]) - &(_matrix[0][0][0][0]));
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0) || (t >= _attr._t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::VoxelToIndex: parameter out of range\n";
    return 0;
  } else {
    return (&(_matrix[t][z][y][x]) - &(_matrix[0][0][0][0]));
  }
#endif
}

template <class VoxelType> inline VoxelType *irtkGenericImage<VoxelType>::GetPointerToVoxels(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return &(_matrix[t][z][y][x]);
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0) || (t >= _attr._t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::GetPointerToVoxels: parameter out of range\n";
    cout << x << " " << y << " " << z << " " << t << endl;
    return NULL;
  } else {
    return &(_matrix[t][z][y][x]);
  }
#endif
}

template <class VoxelType> inline void *irtkGenericImage<VoxelType>::GetScalarPointer(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return &(_matrix[t][z][y][x]);
#else
  if ((x >= _attr._x) || (x < 0) || (y >= _attr._y) || (y < 0) || (z >= _attr._z) || (z < 0) || (t >= _attr._t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::GetScalarPointer: parameter out of range\n";
    cout << x << " " << y << " " << z << " " << t << endl;
    return NULL;
  } else {
    return &(_matrix[t][z][y][x]);
  }
#endif
}

template <> inline int irtkGenericImage<char>::GetScalarType() const
{
	return IRTK_VOXEL_CHAR;
}

template <> inline int irtkGenericImage<unsigned char>::GetScalarType() const
{
	return IRTK_VOXEL_UNSIGNED_CHAR;
}

template <> inline int irtkGenericImage< unsigned short>::GetScalarType() const
{
	return IRTK_VOXEL_UNSIGNED_SHORT;
}

template <> inline int irtkGenericImage<short>::GetScalarType() const
{
	return IRTK_VOXEL_SHORT;
}

template <> inline int irtkGenericImage<int>::GetScalarType() const
{
	return IRTK_VOXEL_INT;
}

template <> inline int irtkGenericImage<unsigned int>::GetScalarType() const
{
	return IRTK_VOXEL_UNSIGNED_INT;
}

template <> inline int irtkGenericImage<float>::GetScalarType() const
{
	return IRTK_VOXEL_FLOAT;
}

template <> inline int irtkGenericImage<double>::GetScalarType() const
{
	return IRTK_VOXEL_DOUBLE;
}

template <class VoxelType> inline double irtkGenericImage<VoxelType>::GetScalarTypeMin() const
{
	return std::numeric_limits<VoxelType>::min();
}

template <class VoxelType> inline double irtkGenericImage<VoxelType>::GetScalarTypeMax() const
{
	return std::numeric_limits<VoxelType>::max();
}

#endif

