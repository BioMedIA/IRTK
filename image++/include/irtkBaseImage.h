/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKBASEIMAGE_H

#define _IRTKBASEIMAGE_H

/**
 * Abstract generic class for 2D or 3D images
 *
 * This abstract base class provides generic support for 2D and 3D image
 * classes. It provides functions for accessing image dimension and origin
 * as well as voxel dimensions. It also provides functions for conversion
 * between image and world coordinates.
 */

class irtkBaseImage : public irtkObject
{

protected:

  /// Image x-dimension (in voxels)
  int _x;
  /// Image y-dimension (in voxels)
  int _y;
  /// Image z-dimension (in voxels)
  int _z;
  /// Image t-dimension (in voxels)
  int _t;

  /// Voxel x-dimensions(in mm)
  double _dx;
  /// Voxel y-dimensions(in mm)
  double _dy;
  /// Voxel z-dimensions(in mm)
  double _dz;
  /// Voxel t-dimensions (in ms)
  double _dt;

  /// Image origin (spatial)
  irtkPoint _origin;

  /// Image origin (temporal)
  double _torigin;

  /// Direction of x-axis
  double _xaxis[3];

  /// Direction of y-axis
  double _yaxis[3];

  /// Direction of z-axis
  double _zaxis[3];

  /// Transformation matrix from image coordinates to world coordinates
  irtkMatrix _matI2W;

  /// Transformation matrix from world coordinates to image coordinates
  irtkMatrix _matW2I;

public:

  /// Default constructor
  irtkBaseImage();

  /// Constructor
  irtkBaseImage(int, int, int, int = 1);

  /// Constructor with given size and dimensions (one frame)
  irtkBaseImage(int, int, int, double, double, double);

  /// Constructor with given size and dimensions (multiple frames)
  irtkBaseImage(int, int, int, int, double, double, double, double);

  /// Constructor with given size, dimensions and axes (one frame)
  irtkBaseImage(int, int, int, double, double, double, double *, double *, double * = NULL);

  /// Constructor with given size, dimensions and axes (multiple frames)
  irtkBaseImage(int, int, int, int, double, double, double, double, double *, double *, double * = NULL);

  /// Constructor with given size, dimensions, origin and axes (one frame)
  irtkBaseImage(int, int, int, double, double, double, irtkPoint, double *, double *, double * = NULL);

  /// Constructor with given size, dimensions, origin and axes (multiple frames)
  irtkBaseImage(int, int, int, int, double, double, double, double, irtkPoint, double, double *, double *, double * = NULL);

  /// Copy constructor
  irtkBaseImage(const irtkBaseImage &);

  /// Destructor
  ~irtkBaseImage();

  /// Update transformation matrix
  void UpdateMatrix();

  /// Initialize baseimage
  void Initialize(int, int, int, int = 1);

  /// Initialize baseimage with given size and dimensions (one frame)
  void Initialize(int, int, int, double, double, double);

  /// Initialize baseimage with given size and dimensions (multiple frames)
  void Initialize(int, int, int, int, double, double, double, double);

  /// Initialize baseimage with given size, dimensions and axes (one frame)
  void Initialize(int, int, int, double, double, double,
                  const double *, const double *, const double * = NULL);

  /// Initialize baseimage with given size, dimensions and axes (multiple frames)
  void Initialize(int, int, int, int, double, double, double, double,
                  const double *, const double *, const double * = NULL);

  /// Initialize baseimage with given size, dimensions, origin and axes (single frame)
  void Initialize(int, int, int, double, double, double, irtkPoint,
                  const double *, const double *, const double * = NULL);

  /// Initialize baseimage with given size, dimensions, origin and axes (multiple frames)
  void Initialize(int, int, int, int, double, double, double, double, irtkPoint, double,
                  const double *, const double *, const double * = NULL);

  /// Initialize baseimage with given size, dimensions and origin
  void Initialize(const irtkBaseImage &);

  /// Copy operators
  irtkBaseImage& operator= (const irtkBaseImage &);

  /// Comparison Operator == (explicit negation replaces != operator)
  Bool       operator==(const irtkBaseImage &);

  //
  // Access functions for image dimensions
  //

  /// Returns the number of voxels in the x-direction
  int  GetX() const;

  /// Returns the number of voxels in the y-direction
  int  GetY() const;

  /// Returns the number of voxels in the z-direction
  int  GetZ() const;

  /// Returns the number of voxels in the t-direction
  int  GetT() const;

  /// Returns the total number of voxels
  int  GetNumberOfVoxels() const;

  //
  // Access functions for voxel dimensions
  //

  /// Returns the number of voxels in the x-direction
  double GetXSize() const;

  /// Returns the number of voxels in the y-direction
  double GetYSize() const;

  /// Returns the number of voxels in the z-direction
  double GetZSize() const;

  /// Returns the number of voxels in the t-direction
  double GetTSize() const;

  /// Voxel dimensions get access
  void  GetPixelSize(double *, double *, double *) const;

  /// Voxel dimensions get access
  void  GetPixelSize(double *, double *, double *, double *) const;

  /// Voxel dimensions put access
  void  PutPixelSize(double, double, double);

  /// Voxel dimensions put access
  void  PutPixelSize(double, double, double, double);

  /// Image origin get access
  irtkPoint GetOrigin() const;

  /// Image origin get access
  void  GetOrigin(double &, double &, double &) const;

  /// Image origin get access
  void  GetOrigin(double &, double &, double &, double &) const;

  /// Image origin put access
  void  PutOrigin(const irtkPoint &);

  /// Image origin put access
  void  PutOrigin(double, double, double);

  /// Image origin put access
  void  PutOrigin(double, double, double, double);

  /// Put image x- and y-axis and z-axis
  void  PutOrientation(double *, double *, double *);

  /// Get image x- and y-axis and z-axis
  void  GetOrientation(double *, double *, double *) const;

  /// Get orientation of axis relative to patient
  void  Orientation(int &, int &, int &) const;

  /// Image to world coordinate conversion with a given point
  void ImageToWorld(irtkPoint &) const;

  /// World to image coordinate conversion with a given point
  void WorldToImage(irtkPoint &) const;

  /// Image to world coordinate conversion with three doubles
  void ImageToWorld(double &, double &, double &) const;

  /// World to image coordinate conversion with three doubles
  void WorldToImage(double &, double &, double &) const;

  /// Return transformation matrix for image to world coordinates
  irtkMatrix GetImageToWorldMatrix() const;

  /// Return transformation matrix for world to image coordinates
  irtkMatrix GetWorldToImageMatrix() const;

  /// Image to time coordinate conversion
  double ImageToTime(double) const;

  /// Time to image coordinate conversion
  double TimeToImage(double) const;

  /// Returns true if point is within the field of view of image
  int IsInFOV(double, double, double);

  /// Print function
  void Print();

  /// Returns the name of the image class
  virtual const char *NameOfClass() = 0;

};

inline int irtkBaseImage::GetX(void) const
{
  return _x;
}

inline int irtkBaseImage::GetY(void) const
{
  return _y;
}

inline int irtkBaseImage::GetZ(void) const
{
  return _z;
}

inline int irtkBaseImage::GetT(void) const
{
  return _t;
}

inline int irtkBaseImage::GetNumberOfVoxels(void) const
{
  return _x*_y*_z*_t;
}

inline double irtkBaseImage::GetXSize(void) const
{
  return _dx;
}

inline double irtkBaseImage::GetYSize(void) const
{
  return _dy;
}

inline double irtkBaseImage::GetZSize(void) const
{
  return _dz;
}

inline double irtkBaseImage::GetTSize(void) const
{
  return _dt;
}

inline void irtkBaseImage::PutPixelSize(double dx, double dy, double dz)
{
  _dx = dx;
  _dy = dy;
  _dz = dz;

  // Update transformation matrix
  this->UpdateMatrix();
}

inline void irtkBaseImage::PutPixelSize(double dx, double dy, double dz, double dt)
{
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = dt;

  // Update transformation matrix
  this->UpdateMatrix();
}

inline void irtkBaseImage::GetPixelSize(double *dx, double *dy, double *dz) const
{
  *dx = _dx;
  *dy = _dy;
  *dz = _dz;
}

inline void irtkBaseImage::GetPixelSize(double *dx, double *dy, double *dz, double *dt) const
{
  *dx = _dx;
  *dy = _dy;
  *dz = _dz;
  *dt = _dt;
}

inline void irtkBaseImage::PutOrigin(const irtkPoint &p)
{
  _origin = p;

  // Update transformation matrix
  this->UpdateMatrix();
}

inline void irtkBaseImage::PutOrigin(double x, double y, double z)
{
  _origin._x = x;
  _origin._y = y;
  _origin._z = z;

  // Update transformation matrix
  this->UpdateMatrix();
}

inline void irtkBaseImage::PutOrigin(double x, double y, double z, double t)
{
  _origin._x = x;
  _origin._y = y;
  _origin._z = z;

  // Update transformation matrix
  this->UpdateMatrix();

  // Calculate origin
  _torigin = t;
}

inline void irtkBaseImage::GetOrigin(double &x, double &y, double &z) const
{
  x = _origin._x;
  y = _origin._y;
  z = _origin._z;
}

inline void irtkBaseImage::GetOrigin(double &x, double &y, double &z, double &t) const
{
  x = _origin._x;
  y = _origin._y;
  z = _origin._z;
  t = _torigin;
}

inline irtkPoint irtkBaseImage::GetOrigin() const
{
  return _origin;
}

inline void irtkBaseImage::PutOrientation(double *xaxis, double *yaxis, double *zaxis)
{
  _xaxis[0] = xaxis[0];
  _xaxis[1] = xaxis[1];
  _xaxis[2] = xaxis[2];
  _yaxis[0] = yaxis[0];
  _yaxis[1] = yaxis[1];
  _yaxis[2] = yaxis[2];
  if (zaxis == NULL) {
    _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
    _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
    _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];
  } else {
    _zaxis[0] = zaxis[0];
    _zaxis[1] = zaxis[1];
    _zaxis[2] = zaxis[2];
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

inline void irtkBaseImage::GetOrientation(double *xaxis, double *yaxis, double *zaxis) const
{
  xaxis[0] = _xaxis[0];
  xaxis[1] = _xaxis[1];
  xaxis[2] = _xaxis[2];
  yaxis[0] = _yaxis[0];
  yaxis[1] = _yaxis[1];
  yaxis[2] = _yaxis[2];
  zaxis[0] = _zaxis[0];
  zaxis[1] = _zaxis[1];
  zaxis[2] = _zaxis[2];
}

inline void irtkBaseImage::ImageToWorld(double &x, double &y, double &z) const
{
  double a, b, c;

  // Pre-multiply point with transformation matrix
  a = _matI2W(0, 0)*x+_matI2W(0, 1)*y+_matI2W(0, 2)*z+_matI2W(0, 3);
  b = _matI2W(1, 0)*x+_matI2W(1, 1)*y+_matI2W(1, 2)*z+_matI2W(1, 3);
  c = _matI2W(2, 0)*x+_matI2W(2, 1)*y+_matI2W(2, 2)*z+_matI2W(2, 3);

  // Copy result back
  x = a;
  y = b;
  z = c;
}

inline void irtkBaseImage::ImageToWorld(irtkPoint &p) const
{
  double a, b, c;

  // Pre-multiply point with transformation matrix
  a = _matI2W(0, 0)*p._x+_matI2W(0, 1)*p._y+_matI2W(0, 2)*p._z+_matI2W(0, 3);
  b = _matI2W(1, 0)*p._x+_matI2W(1, 1)*p._y+_matI2W(1, 2)*p._z+_matI2W(1, 3);
  c = _matI2W(2, 0)*p._x+_matI2W(2, 1)*p._y+_matI2W(2, 2)*p._z+_matI2W(2, 3);

  // Copy result back
  p._x = a;
  p._y = b;
  p._z = c;
}

inline void irtkBaseImage::WorldToImage(double &x, double &y, double &z) const
{
  double a, b, c;

  // Pre-multiply point with transformation matrix
  a = _matW2I(0, 0)*x+_matW2I(0, 1)*y+_matW2I(0, 2)*z+_matW2I(0, 3);
  b = _matW2I(1, 0)*x+_matW2I(1, 1)*y+_matW2I(1, 2)*z+_matW2I(1, 3);
  c = _matW2I(2, 0)*x+_matW2I(2, 1)*y+_matW2I(2, 2)*z+_matW2I(2, 3);

  // Copy result back
  x = a;
  y = b;
  z = c;
}

inline void irtkBaseImage::WorldToImage(irtkPoint &p) const
{
  double a, b, c;

  // Pre-multiply point with transformation matrix
  a = _matW2I(0, 0)*p._x+_matW2I(0, 1)*p._y+_matW2I(0, 2)*p._z+_matW2I(0, 3);
  b = _matW2I(1, 0)*p._x+_matW2I(1, 1)*p._y+_matW2I(1, 2)*p._z+_matW2I(1, 3);
  c = _matW2I(2, 0)*p._x+_matW2I(2, 1)*p._y+_matW2I(2, 2)*p._z+_matW2I(2, 3);

  // Copy result back
  p._x = a;
  p._y = b;
  p._z = c;
}

inline irtkMatrix irtkBaseImage::GetImageToWorldMatrix() const
{
  return _matI2W;
}

inline irtkMatrix irtkBaseImage::GetWorldToImageMatrix() const
{
  return _matW2I;
}

inline double irtkBaseImage::ImageToTime(double t) const
{
  return _torigin+t*_dt;
}

inline double irtkBaseImage::TimeToImage(double t) const
{
  return (t-_torigin)/_dt;
}

inline int irtkBaseImage::IsInFOV(double x, double y, double z)
{
  this->WorldToImage(x, y, z);
  if ((x < -0.5) || (x >= _x-0.5) ||
      (y < -0.5) || (y >= _y-0.5) ||
      (z < -0.5) || (z >= _z-0.5)) {
    return False;
  }
  return True;
}

#endif
