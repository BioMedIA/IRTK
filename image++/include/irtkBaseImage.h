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

#include <irtkImageAttributes.h>

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

	/// Attributes of the imaging geometry;
	irtkImageAttributes _attr;

  /// Transformation matrix from image coordinates to world coordinates
  irtkMatrix _matI2W;

  /// Transformation matrix from world coordinates to image coordinates
  irtkMatrix _matW2I;

  /// Constructor for empty image
  irtkBaseImage();

  /// Constructor for image with certain image attributes
  irtkBaseImage(const irtkImageAttributes &);

  /// Copy constructor
  irtkBaseImage(const irtkBaseImage &);

  /// Update baseimage
  void Update(const irtkImageAttributes &);

  /// Update transformation matrix
  void UpdateMatrix();

public:

  /// Destructor
  virtual ~irtkBaseImage();

  //
  // Access functions for image dimensions
  //

  /// Returns the number of voxels in the x-direction
  virtual int  GetX() const;

  /// Returns the number of voxels in the y-direction
  virtual int  GetY() const;

  /// Returns the number of voxels in the z-direction
  virtual int  GetZ() const;

  /// Returns the number of voxels in the t-direction
  virtual int  GetT() const;

  /// Returns the total number of voxels
  virtual int  GetNumberOfVoxels() const;

  /// Gets the image attributes
  virtual irtkImageAttributes GetImageAttributes() const;
  
  /// Initialize image from attributes
  virtual void Initialize(const irtkImageAttributes &) = 0;

  //
  // Access functions for voxel dimensions
  //

  /// Returns the number of voxels in the x-direction
  virtual double GetXSize() const;

  /// Returns the number of voxels in the y-direction
  virtual double GetYSize() const;

  /// Returns the number of voxels in the z-direction
  virtual double GetZSize() const;

  /// Returns the number of voxels in the t-direction
  virtual double GetTSize() const;

  /// Voxel dimensions get access
  virtual void  GetPixelSize(double *, double *, double *) const;

  /// Voxel dimensions get access
  virtual void  GetPixelSize(double *, double *, double *, double *) const;

  /// Voxel dimensions put access
  virtual void  PutPixelSize(double, double, double);

  /// Voxel dimensions put access
  virtual void  PutPixelSize(double, double, double, double);

  /// Image origin get access
  virtual irtkPoint GetOrigin() const;

  /// Image origin get access
  virtual void  GetOrigin(double &, double &, double &) const;

  /// Image origin get access
  virtual void  GetOrigin(double &, double &, double &, double &) const;

  /// Image origin put access
  virtual void  PutOrigin(const irtkPoint &);

  /// Image origin put access
  virtual void  PutOrigin(double, double, double);

  /// Image origin put access
  virtual void  PutOrigin(double, double, double, double);

  /// Put image x- and y-axis and z-axis
  virtual void  PutOrientation(double *, double *, double *);

  /// Get image x- and y-axis and z-axis
  virtual void  GetOrientation(double *, double *, double *) const;

  /// Get orientation of axis relative to patient
  virtual void  Orientation(int &, int &, int &) const;

  /// Image to world coordinate conversion with a given point
  virtual void ImageToWorld(irtkPoint &) const;

  /// World to image coordinate conversion with a given point
  virtual void WorldToImage(irtkPoint &) const;

  /// Image to world coordinate conversion with three doubles
  virtual void ImageToWorld(double &, double &, double &) const;

  /// World to image coordinate conversion with three doubles
  virtual void WorldToImage(double &, double &, double &) const;

  /// Return transformation matrix for image to world coordinates
  virtual irtkMatrix GetImageToWorldMatrix() const;

  /// Return transformation matrix for world to image coordinates
  virtual irtkMatrix GetWorldToImageMatrix() const;

  /// Image to time coordinate conversion
  virtual double ImageToTime(double) const;

  /// Time to image coordinate conversion
  virtual double TimeToImage(double) const;

  /// Returns true if point is within the field of view of image
  virtual bool IsInFOV(double, double, double);

  /// boolean operation for empty
  virtual bool IsEmpty() const;

  /// Minimum and maximum pixel values get accessor
  virtual void GetMinMaxAsDouble(double *, double *) const;

  /// Minimum and maximum pixel values put accessor
  virtual void PutMinMaxAsDouble(double, double);

  /// Function for pixel get access as double
  virtual double GetAsDouble(int, int, int, int = 0) const = 0;

  /// Function for pixel put access
  virtual void   PutAsDouble(int, int, int, double) = 0;

  /// Function for pixel put access
  virtual void   PutAsDouble(int, int, int, int, double) = 0;

  /// Returns the name of the image class
  virtual const char *NameOfClass() = 0;

  /// Function for pixel access via pointers
  virtual void *GetScalarPointer(int = 0, int = 0, int = 0, int = 0) const = 0;

  /// Function which returns pixel scalar type
  virtual int GetScalarType() const = 0;

  /// Function which returns the minimum value the pixel can hold without overflowing
  virtual double GetScalarTypeMin() const = 0;

  /// Function which returns the minimum value the pixel can hold without overflowing
  virtual double GetScalarTypeMax() const = 0;

  /// Reflect image around x
  virtual void ReflectX() = 0;

  /// Reflect image around y
  virtual void ReflectY() = 0;

  /// Reflect image around z
  virtual void ReflectZ() = 0;

  /// Flip x and y axis
  virtual void FlipXY() = 0;

  /// Flip x and z axis
  virtual void FlipXZ() = 0;

  /// Flip y and z axis
  virtual void FlipYZ() = 0;

  /// Flip x and t axis
  virtual void FlipXT() = 0;

  /// Flip y and t axis
  virtual void FlipYT() = 0;

  /// Flip z and t axis
  virtual void FlipZT() = 0;

  /// Read file and construct image
  static irtkBaseImage *New(const char *);

  /// Read file and construct image
  static irtkBaseImage *New(const irtkBaseImage *);

  /// Write file
  virtual void Write(const char *) = 0;

  /// Print function
  virtual void Print();

};

inline int irtkBaseImage::GetX() const
{
  return _attr._x;
}

inline int irtkBaseImage::GetY() const
{
  return _attr._y;
}

inline int irtkBaseImage::GetZ() const
{
  return _attr._z;
}

inline int irtkBaseImage::GetT() const
{
  return _attr._t;
}

inline irtkImageAttributes irtkBaseImage::GetImageAttributes() const
{
  return _attr;
}

inline int irtkBaseImage::GetNumberOfVoxels(void) const
{
  return _attr._x*_attr._y*_attr._z*_attr._t;
}

inline double irtkBaseImage::GetXSize(void) const
{
  return _attr._dx;
}

inline double irtkBaseImage::GetYSize(void) const
{
  return _attr._dy;
}

inline double irtkBaseImage::GetZSize(void) const
{
  return _attr._dz;
}

inline double irtkBaseImage::GetTSize(void) const
{
  return _attr._dt;
}

inline void irtkBaseImage::PutPixelSize(double dx, double dy, double dz)
{
	_attr._dx = dx;
	_attr._dy = dy;
	_attr._dz = dz;

  // Update transformation matrix
  this->UpdateMatrix();
}

inline void irtkBaseImage::PutPixelSize(double dx, double dy, double dz, double dt)
{
	_attr._dx = dx;
	_attr._dy = dy;
	_attr._dz = dz;
	_attr._dt = dt;

  // Update transformation matrix
  this->UpdateMatrix();
}

inline void irtkBaseImage::GetPixelSize(double *dx, double *dy, double *dz) const
{
  *dx = _attr._dx;
  *dy = _attr._dy;
  *dz = _attr._dz;
}

inline void irtkBaseImage::GetPixelSize(double *dx, double *dy, double *dz, double *dt) const
{
  *dx = _attr._dx;
  *dy = _attr._dy;
  *dz = _attr._dz;
  *dt = _attr._dt;
}

inline void irtkBaseImage::PutOrigin(const irtkPoint &p)
{
	_attr._xorigin = p._x;
	_attr._yorigin = p._y;
	_attr._zorigin = p._z;

	// Update transformation matrix
	this->UpdateMatrix();
}

inline void irtkBaseImage::PutOrigin(double x, double y, double z)
{
	_attr._xorigin = x;
	_attr._yorigin = y;
	_attr._zorigin = z;

	// Update transformation matrix
	this->UpdateMatrix();
}

inline void irtkBaseImage::PutOrigin(double x, double y, double z, double t)
{
	_attr._xorigin = x;
	_attr._yorigin = y;
	_attr._zorigin = z;
	// Calculate origin
	_attr._torigin = t;

	// Update transformation matrix
	this->UpdateMatrix();
}

inline void irtkBaseImage::GetOrigin(double &x, double &y, double &z) const
{
  x = _attr._xorigin;
  y = _attr._yorigin;
  z = _attr._zorigin;
}

inline void irtkBaseImage::GetOrigin(double &x, double &y, double &z, double &t) const
{
  x = _attr._xorigin;
  y = _attr._yorigin;
  z = _attr._zorigin;
  t = _attr._torigin;
}

inline irtkPoint irtkBaseImage::GetOrigin() const
{
  return irtkPoint(_attr._xorigin, _attr._yorigin, _attr._zorigin);
}

inline void irtkBaseImage::PutOrientation(double *xaxis, double *yaxis, double *zaxis)
{
	_attr._xaxis[0] = xaxis[0];
	_attr._xaxis[1] = xaxis[1];
	_attr._xaxis[2] = xaxis[2];
	_attr._yaxis[0] = yaxis[0];
	_attr._yaxis[1] = yaxis[1];
	_attr._yaxis[2] = yaxis[2];
  if (zaxis == NULL) {
  	_attr._zaxis[0] = _attr._xaxis[1]*_attr._yaxis[2] - _attr._xaxis[2]*_attr._yaxis[1];
  	_attr._zaxis[1] = _attr._xaxis[2]*_attr._yaxis[0] - _attr._xaxis[0]*_attr._yaxis[2];
  	_attr._zaxis[2] = _attr._xaxis[0]*_attr._yaxis[1] - _attr._xaxis[1]*_attr._yaxis[0];
  } else {
  	_attr._zaxis[0] = zaxis[0];
  	_attr._zaxis[1] = zaxis[1];
  	_attr._zaxis[2] = zaxis[2];
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

inline void irtkBaseImage::GetOrientation(double *xaxis, double *yaxis, double *zaxis) const
{
  xaxis[0] = _attr._xaxis[0];
  xaxis[1] = _attr._xaxis[1];
  xaxis[2] = _attr._xaxis[2];
  yaxis[0] = _attr._yaxis[0];
  yaxis[1] = _attr._yaxis[1];
  yaxis[2] = _attr._yaxis[2];
  zaxis[0] = _attr._zaxis[0];
  zaxis[1] = _attr._zaxis[1];
  zaxis[2] = _attr._zaxis[2];
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
  return _attr._torigin+t*_attr._dt;
}

inline double irtkBaseImage::TimeToImage(double t) const
{
  return (t-_attr._torigin)/_attr._dt;
}

inline bool irtkBaseImage::IsInFOV(double x, double y, double z)
{
  this->WorldToImage(x, y, z);
  if ((x < -0.5) || (x >= _attr._x-0.5) ||
      (y < -0.5) || (y >= _attr._y-0.5) ||
      (z < -0.5) || (z >= _attr._z-0.5)) {
    return false;
  }
  return true;
}

inline bool irtkBaseImage::IsEmpty() const
{
  return ((_attr._x < 1) || (_attr._y < 1) || (_attr._z < 1) || (_attr._t < 1));
}

typedef class irtkBaseImage irtkImage;

#endif
