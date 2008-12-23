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
  irtkGenericImage(void);

  /// Constructor from image file
  irtkGenericImage(char *);

  /// Constructor for given image dimensions
  irtkGenericImage(int, int, int, int = 1);

  /// Constructor for given image dimensions and voxel dimensions (single frame)
  irtkGenericImage(int, int, int, double, double, double);

  /// Constructor for given image dimensions and voxel dimensions (multiple frames)
  irtkGenericImage(int, int, int, int, double, double, double, double);

  /// Constructor for given image dimensions, voxel dimensions and orientation (single frame)
  irtkGenericImage(int, int, int, double, double, double,
                   const double *, const double *, const double * = NULL);

  /// Constructor for given image dimensions, voxel dimensions and orientation (multiple frames)
  irtkGenericImage(int, int, int, int, double, double, double, double,
                   const double *, const double *, const double * = NULL);

  /// Constructor for given image dimensions, voxel dimensions, origin (single frame)
  /// and orientation
  irtkGenericImage(int, int, int, double, double, double, irtkPoint,
                   const double *, const double *, const double * = NULL);

  /// Constructor for given image dimensions, voxel dimensions, origin (multiple frames)
  /// and orientation
  irtkGenericImage(int, int, int, int, double, double, double, double, irtkPoint, double,
                   const double *, const double *, const double * = NULL);

  /// Copy constructor for image of type unsigned char
  irtkGenericImage(const irtkGenericImage<irtkBytePixel> &);

  /// Copy constructor for image of type short
  irtkGenericImage(const irtkGenericImage<irtkGreyPixel> &);

  /// Copy constructor for image of type float
  irtkGenericImage(const irtkGenericImage<irtkRealPixel> &);

  /// Copy constructor for image of type RGB pixel
  irtkGenericImage(const irtkGenericImage<irtkRGBPixel> &);

  /// Copy constructor for image of type irtkVector3D<char>.
  irtkGenericImage(const irtkGenericImage<irtkVector3D<char> > &);

  /// Copy constructor for image of type irtkVector3D<short>.
  irtkGenericImage(const irtkGenericImage<irtkVector3D<short> > &);

  /// Copy constructor for image of type irtkVector3D<float>.
  irtkGenericImage(const irtkGenericImage<irtkVector3D<float> > &);

  /// Copy constructor for image of type irtkVector3D<double>.
  irtkGenericImage(const irtkGenericImage<irtkVector3D<double> > &);

  /// Destructor
  ~irtkGenericImage(void);

  /// Initialize an image
  void Initialize(int, int, int, int = 1);

  /// Initialize baseimage with given size and dimensions (single frame)
  void Initialize(int, int, int, double, double, double);

  /// Initialize baseimage with given size and dimensions (multiple frames)
  void Initialize(int, int, int, int, double, double, double, double);

  /// Initialize baseimage with given size, dimensions and axes (single frame)
  void Initialize(int, int, int, double, double, double,
                  const double *, const double *, const double * = NULL);

  /// Initialize baseimage with given size, dimensions and axes (multiple frames)
  void Initialize(int, int, int, int, double, double, double, double,
                  const double *, const double *, const double * = NULL);

  /// Initialize baseimage with given size, dimensions, origin and axes (single frame)
  void Initialize(int, int, int, double, double, double, irtkPoint,
                  const double *, const double *, const double * = NULL);

  /// Initialize baseimage with given size, dimensions, origin and axes (multiple frames)
  void Initialize(int, int, int, int, double, double, double, double, irtkPoint, double, const double *, const double *, const double * = NULL);

  /// Initialize an image
  void Initialize(const irtkBaseImage &);

  /// Read image from file
  void Read (const char *);

  /// Write image to file
  void Write(const char *);

  /// Minimum and maximum pixel values get accessor
  void GetMinMax(VoxelType *, VoxelType *) const;

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

  /// Function for pixel get access as double
  double GetAsDouble(int, int, int, int = 0) const;

  /// Function for pixel put access
  void   Put(int, int, int, VoxelType);

  /// Function for pixel put access
  void   Put(int, int, int, int, VoxelType);

  /// Function for pixel put access
  void   PutAsDouble(int, int, int, double);

  /// Function for pixel put access
  void   PutAsDouble(int, int, int, int, double);

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

  /// Copy operator for image of type unsigned char
  irtkGenericImage& operator= (const irtkGenericImage<irtkBytePixel> &);
  /// Copy operator for image of type short
  irtkGenericImage& operator= (const irtkGenericImage<irtkGreyPixel> &);
  /// Copy operator for image of type float
  irtkGenericImage& operator= (const irtkGenericImage<irtkRealPixel> &);
  /// Copy operator for image of type RGB pixel
  irtkGenericImage& operator= (const irtkGenericImage<irtkRGBPixel> &);
  /// Assignment operator for image of type irtkVector3D<char>.
  irtkGenericImage& operator= (const irtkGenericImage<irtkVector3D<char> > &);
  /// Assignment operator for image of type irtkVector3D<short>.
  irtkGenericImage& operator= (const irtkGenericImage<irtkVector3D<short> > &);
  /// Assignment operator for image of type irtkVector3D<float>.
  irtkGenericImage& operator= (const irtkGenericImage<irtkVector3D<float> > &);
  /// Assignment operator for image of type irtkVector3D<double>.
  irtkGenericImage& operator= (const irtkGenericImage<irtkVector3D<double> > &);

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
  Bool operator==(const irtkGenericImage &);
  /// Comparison operator != (if _HAS_STL is defined, negate == operator)
  ///  Bool operator!=(const irtkGenericImage &);

  /// Boolean operation for empty
  Bool IsEmpty() const;

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
  void FlipXY();
  /// Flip x and z axis
  void FlipXZ();
  /// Flip y and z axis
  void FlipYZ();
  /// Flip x and t axis
  void FlipXT();
  /// Flip y and t axis
  void FlipYT();
  /// Flip z and t axis
  void FlipZT();

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

  /// Returns the name of the image class
  const char *NameOfClass();

};

template <class VoxelType> inline void irtkGenericImage<VoxelType>::Put(int x, int y, int z, VoxelType val)
{
#ifdef NO_BOUNDS
  _matrix[0][z][y][x] = val;
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (_t != 0)) {
    cout << "irtkGenericImage<VoxelType>::Put: parameter out of range\n";
  } else {
    _matrix[0][z][y][x] = val;
  }
#endif
}

template <class VoxelType> inline void irtkGenericImage<VoxelType>::Put(int x, int y, int z, int t, VoxelType val)
{
#ifdef NO_BOUNDS
  _matrix[t][z][y][x] = val;
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<VoxelType>::Put: parameter out of range\n";
  } else {
    _matrix[t][z][y][x] = val;
  }
#endif
}

template <> inline void irtkGenericImage<irtkBytePixel>::PutAsDouble(int x, int y, int z, double val)
{
  if (val > MAX_BYTE) val = MAX_BYTE;
  if (val < MIN_BYTE) val = MIN_BYTE;
#ifdef NO_BOUNDS
  _matrix[0][z][y][x] = round(val);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (_t != 0)) {
    cout << "irtkGenericImage<Type>::PutAsDouble: parameter out of range\n";
  } else {
    _matrix[0][z][y][x] = round(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkBytePixel>::PutAsDouble(int x, int y, int z, int t, double val)
{
  if (val > MAX_BYTE) val = MAX_BYTE;
  if (val < MIN_BYTE) val = MIN_BYTE;
#ifdef NO_BOUNDS
  _matrix[t][z][y][x] = round(val);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::PutAsDouble: parameter out of range\n";
  } else {
    _matrix[t][z][y][x] = round(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkGreyPixel>::PutAsDouble(int x, int y, int z, double val)
{
  if (val > MAX_GREY) val = MAX_GREY;
  if (val < MIN_GREY) val = MIN_GREY;
#ifdef NO_BOUNDS
  _matrix[0][z][y][x] = round(val);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (_t != 0)) {
    cout << "irtkGenericImage<Type>::PutAsDouble: parameter out of range\n";
  } else {
    _matrix[0][z][y][x] = round(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkGreyPixel>::PutAsDouble(int x, int y, int z, int t, double val)
{
  if (val > MAX_GREY) val = MAX_GREY;
  if (val < MIN_GREY) val = MIN_GREY;
#ifdef NO_BOUNDS
  _matrix[t][z][y][x] = round(val);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::PutAsDouble: parameter out of range\n";
  } else {
    _matrix[t][z][y][x] = round(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkRealPixel>::PutAsDouble(int x, int y, int z, double val)
{
#ifdef NO_BOUNDS
  _matrix[0][z][y][x] = static_cast<irtkRealPixel>(val);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (_t != 0)) {
    cout << "irtkGenericImage<Type>::PutAsDouble: parameter out of range\n";
  } else {
    _matrix[0][z][y][x] = static_cast<irtkRealPixel>(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkRealPixel>::PutAsDouble(int x, int y, int z, int t, double val)
{
#ifdef NO_BOUNDS
  _matrix[t][z][y][x] = static_cast<irtkRealPixel>(val);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::PutAsDouble: parameter out of range\n";
  } else {
    _matrix[t][z][y][x] = static_cast<irtkRealPixel>(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkVector3D<char> >::PutAsDouble(int x, int y, int z, double val)
{
#ifdef NO_BOUNDS
  val = val/sqrt(3.0);
  _matrix[0][z][y][x]._x = static_cast<char>(val);
  _matrix[0][z][y][x]._y = static_cast<char>(val);
  _matrix[0][z][y][x]._z = static_cast<char>(val);
#else
  if (x >= _x || x < 0 || y >= _y || y < 0 || z >= _z || z < 0 || _t != 0) {
    cout << "irtkGenericImage<irtkVector3D<char> >::PutAsDouble: parameter out of range\n";
  } else {
    val = val/sqrt(3.0);
    _matrix[0][z][y][x]._x = static_cast<char>(val);
    _matrix[0][z][y][x]._y = static_cast<char>(val);
    _matrix[0][z][y][x]._z = static_cast<char>(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkVector3D<char> >::PutAsDouble(int x, int y, int z, int t, double val)
{
#ifdef NO_BOUNDS
  val = val/sqrt(3.0);
  _matrix[t][z][y][x]._x = static_cast<char>(val);
  _matrix[t][z][y][x]._y = static_cast<char>(val);
  _matrix[t][z][y][x]._z = static_cast<char>(val);
#else
  if (x >= _x || x < 0 || y >= _y || y < 0 || z >= _z || z < 0 || t >= _t || t < 0) {
    cout << "irtkGenericImage<irtkVector3D<char> >::PutAsDouble: parameter out of range\n";
  } else {
    val = val/sqrt(3.0);
    _matrix[t][z][y][x]._x = static_cast<char>(val);
    _matrix[t][z][y][x]._y = static_cast<char>(val);
    _matrix[t][z][y][x]._z = static_cast<char>(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkVector3D<short> >::PutAsDouble(int x, int y, int z, double val)
{
#ifdef NO_BOUNDS
  val = val/sqrt(3.0);
  _matrix[0][z][y][x]._x = static_cast<short>(val);
  _matrix[0][z][y][x]._y = static_cast<short>(val);
  _matrix[0][z][y][x]._z = static_cast<short>(val);
#else
  if (x >= _x || x < 0 || y >= _y || y < 0 || z >= _z || z < 0 || _t != 0) {
    cout << "irtkGenericImage<irtkVector3D<short> >::PutAsDouble: parameter out of range\n";
  } else {
    val = val/sqrt(3.0);
    _matrix[0][z][y][x]._x = static_cast<short>(val);
    _matrix[0][z][y][x]._y = static_cast<short>(val);
    _matrix[0][z][y][x]._z = static_cast<short>(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkVector3D<short> >::PutAsDouble(int x, int y, int z, int t, double val)
{
#ifdef NO_BOUNDS
  val = val/sqrt(3.0);
  _matrix[t][z][y][x]._x = static_cast<short>(val);
  _matrix[t][z][y][x]._y = static_cast<short>(val);
  _matrix[t][z][y][x]._z = static_cast<short>(val);
#else
  if (x >= _x || x < 0 || y >= _y || y < 0 || z >= _z || z < 0 || t >= _t || t < 0) {
    cout << "irtkGenericImage<irtkVector3D<short> >::PutAsDouble: parameter out of range\n";
  } else {
    val = val/sqrt(3.0);
    _matrix[t][z][y][x]._x = static_cast<short>(val);
    _matrix[t][z][y][x]._y = static_cast<short>(val);
    _matrix[t][z][y][x]._z = static_cast<short>(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkVector3D<float> >::PutAsDouble(int x, int y, int z, double val)
{
#ifdef NO_BOUNDS
  val = val/sqrt(3.0);
  _matrix[0][z][y][x]._x = static_cast<float>(val);
  _matrix[0][z][y][x]._y = static_cast<float>(val);
  _matrix[0][z][y][x]._z = static_cast<float>(val);
#else
  if (x >= _x || x < 0 || y >= _y || y < 0 || z >= _z || z < 0 || _t != 0) {
    cout << "irtkGenericImage<irtkVector3D<float> >::PutAsDouble: parameter out of range\n";
  } else {
    val = val/sqrt(3.0);
    _matrix[0][z][y][x]._x = static_cast<float>(val);
    _matrix[0][z][y][x]._y = static_cast<float>(val);
    _matrix[0][z][y][x]._z = static_cast<float>(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkVector3D<float> >::PutAsDouble(int x, int y, int z, int t, double val)
{
#ifdef NO_BOUNDS
  val = val/sqrt(3.0);
  _matrix[t][z][y][x]._x = static_cast<float>(val);
  _matrix[t][z][y][x]._y = static_cast<float>(val);
  _matrix[t][z][y][x]._z = static_cast<float>(val);
#else
  if (x >= _x || x < 0 || y >= _y || y < 0 || z >= _z || z < 0 || t >= _t || t < 0) {
    cout << "irtkGenericImage<irtkVector3D<float> >::PutAsDouble: parameter out of range\n";
  } else {
    val = val/sqrt(3.0);
    _matrix[t][z][y][x]._x = static_cast<float>(val);
    _matrix[t][z][y][x]._y = static_cast<float>(val);
    _matrix[t][z][y][x]._z = static_cast<float>(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkVector3D<double> >::PutAsDouble(int x, int y, int z, double val)
{
#ifdef NO_BOUNDS
  val = val/sqrt(3.0);
  _matrix[0][z][y][x]._x = static_cast<double>(val);
  _matrix[0][z][y][x]._y = static_cast<double>(val);
  _matrix[0][z][y][x]._z = static_cast<double>(val);
#else
  if (x >= _x || x < 0 || y >= _y || y < 0 || z >= _z || z < 0 || _t != 0) {
    cout << "irtkGenericImage<irtkVector3D<double> >::PutAsDouble: parameter out of range\n";
  } else {
    val = val/sqrt(3.0);
    _matrix[0][z][y][x]._x = static_cast<double>(val);
    _matrix[0][z][y][x]._y = static_cast<double>(val);
    _matrix[0][z][y][x]._z = static_cast<double>(val);
  }
#endif
}

template <> inline void irtkGenericImage<irtkVector3D<double> >::PutAsDouble(int x, int y, int z, int t, double val)
{
#ifdef NO_BOUNDS
  val = val/sqrt(3.0);
  _matrix[t][z][y][x]._x = static_cast<double>(val);
  _matrix[t][z][y][x]._y = static_cast<double>(val);
  _matrix[t][z][y][x]._z = static_cast<double>(val);
#else
  if (x >= _x || x < 0 || y >= _y || y < 0 || z >= _z || z < 0 || t >= _t || t < 0) {
    cout << "irtkGenericImage<irtkVector3D<double> >::PutAsDouble: parameter out of range\n";
  } else {
    val = val/sqrt(3.0);
    _matrix[t][z][y][x]._x = static_cast<double>(val);
    _matrix[t][z][y][x]._y = static_cast<double>(val);
    _matrix[t][z][y][x]._z = static_cast<double>(val);
  }
#endif
}

template <class VoxelType> inline VoxelType irtkGenericImage<VoxelType>::Get(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return (_matrix[t][z][y][x]);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (z < t)) {
    cout << "irtkGenericImage<Type>::Get: parameter out of range\n";
    return 0;
  } else {
    return(_matrix[t][z][y][x]);
  }
#endif
}

template <> inline double irtkGenericImage<irtkBytePixel>::GetAsDouble(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return (_matrix[t][z][y][x]);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::GetAsDouble: parameter out of range\n";
    return 0;
  } else {
    return (_matrix[t][z][y][x]);
  }
#endif

}

template <> inline double irtkGenericImage<irtkGreyPixel>::GetAsDouble(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return (_matrix[t][z][y][x]);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::GetAsDouble: parameter out of range\n";
    return 0;
  } else {
    return (_matrix[t][z][y][x]);
  }
#endif
}

template <> inline double irtkGenericImage<irtkRealPixel>::GetAsDouble(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return (_matrix[t][z][y][x]);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::GetAsDouble: parameter out of range\n";
    return 0;
  } else {
    return (_matrix[t][z][y][x]);
  }
#endif
}

template <> inline double irtkGenericImage<irtkVector3D<char> >::GetAsDouble(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return sqrt(static_cast<double>(_matrix[t][z][y][x]._x*_matrix[t][z][y][x]._x +
                                  _matrix[t][z][y][x]._y*_matrix[t][z][y][x]._y +
                                  _matrix[t][z][y][x]._z*_matrix[t][z][y][x]._z));
#else
  if ((x >= x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<irtkVector3D<char> >::GetAsDouble: parameter out of range\n";
    return 0;
  } else {
    return sqrt(static_cast<double>(_matrix[t][z][y][x]._x*_matrix[t][z][y][x]._x +
                                    _matrix[t][z][y][x]._y*_matrix[t][z][y][x]._y +
                                    _matrix[t][z][y][x]._z*_matrix[t][z][y][x]._z));
  }
#endif
}

template <> inline double irtkGenericImage<irtkVector3D<short> >::GetAsDouble(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return sqrt(static_cast<double>(_matrix[t][z][y][x]._x*_matrix[t][z][y][x]._x +
                                  _matrix[t][z][y][x]._y*_matrix[t][z][y][x]._y +
                                  _matrix[t][z][y][x]._z*_matrix[t][z][y][x]._z));
#else
  if ((x >= x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<irtkVector3D<short> >::GetAsDouble: parameter out of range\n";
    return 0;
  } else {
    return sqrt(static_cast<double>(_matrix[t][z][y][x]._x*_matrix[t][z][y][x]._x +
                                    _matrix[t][z][y][x]._y*_matrix[t][z][y][x]._y +
                                    _matrix[t][z][y][x]._z*_matrix[t][z][y][x]._z));
  }
#endif
}

template <> inline double irtkGenericImage<irtkVector3D<float> >::GetAsDouble(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return sqrt(static_cast<double>(_matrix[t][z][y][x]._x*_matrix[t][z][y][x]._x +
                                  _matrix[t][z][y][x]._y*_matrix[t][z][y][x]._y +
                                  _matrix[t][z][y][x]._z*_matrix[t][z][y][x]._z));
#else
  if ((x >= x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<irtkVector3D<float> >::GetAsDouble: parameter out of range\n";
    return 0;
  } else {
    return sqrt(static_cast<double>(_matrix[t][z][y][x]._x*_matrix[t][z][y][x]._x +
                                    _matrix[t][z][y][x]._y*_matrix[t][z][y][x]._y +
                                    _matrix[t][z][y][x]._z*_matrix[t][z][y][x]._z));
  }
#endif
}

template <> inline double irtkGenericImage<irtkVector3D<double> >::GetAsDouble(int x, int y, int z, int t) const
{
#ifdef NO_BOUNDS
  return sqrt(static_cast<double>(_matrix[t][z][y][x]._x*_matrix[t][z][y][x]._x +
                                  _matrix[t][z][y][x]._y*_matrix[t][z][y][x]._y +
                                  _matrix[t][z][y][x]._z*_matrix[t][z][y][x]._z));
#else
  if ((x >= x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<irtkVector3D<double> >::GetAsDouble: parameter out of range\n";
    return 0;
  } else {
    return sqrt(static_cast<double>(_matrix[t][z][y][x]._x*_matrix[t][z][y][x]._x +
                                    _matrix[t][z][y][x]._y*_matrix[t][z][y][x]._y +
                                    _matrix[t][z][y][x]._z*_matrix[t][z][y][x]._z));
  }
#endif
}

template <class VoxelType> inline VoxelType& irtkGenericImage<VoxelType>::operator()(int x, int y, int z, int t)
{
#ifdef NO_BOUNDS
  return (_matrix[t][z][y][x]);
#else
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
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
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
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
  if ((x >= _x) || (x < 0) || (y >= _y) || (y < 0) || (z >= _z) || (z < 0) || (t >= _t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::GetPointerToVoxels: parameter out of range\n";
    cout << x << " " << y << " " << z << " " << t << endl;
    return NULL;
  } else {
    return &(_matrix[t][z][y][x]);
  }
#endif
}

template <class VoxelType> inline Bool irtkGenericImage<VoxelType>::IsEmpty() const
{
  return ((_x < 1) || (_y < 1) || (_z < 1) || (_t < 1));
}

#endif

