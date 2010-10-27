/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKFREEFORMTRANSFORMATION3D_H

#define _IRTKFREEFORMTRANSFORMATION3D_H

#include <irtkGeometry.h>

#define FFDLOOKUPTABLESIZE 1000

/**
 * Class for 3D free form transformations
 *
 * This class implements a 3D free form transformation
 *
 */

class irtkFreeFormTransformation3D : public irtkFreeFormTransformation
{

protected:

  /// Number of control points in x
  int _x;

  /// Number of control points in y
  int _y;

  /// Number of control points in z
  int _z;

  /// Spacing of control points in x (in mm)
  double _dx;

  /// Spacing of control points in y (in mm)
  double _dy;

  /// Spacing of control points in z (in mm)
  double _dz;

  /// Direction of x-axis
  double _xaxis[3];

  /// Direction of y-axis
  double _yaxis[3];

  /// Direction of z-axis
  double _zaxis[3];

  /// Origin
  irtkPoint _origin;

  /// Transformation matrix from lattice coordinates to world coordinates
  irtkMatrix _matL2W;

  /// Transformation matrix from world coordinates to lattice coordinates
  irtkMatrix _matW2L;

  /// Displacement in the x-direction at the control points (in mm)
  double ***_xdata;

  /// Displacement in the y-direction at the control points (in mm)
  double ***_ydata;

  /// Displacement in the z-direction at the control points (in mm)
  double ***_zdata;

  /// Allocate memory for control points
  static double ***Allocate  (double ***, int, int, int);

  /// Deallocate memory for control points
  static double ***Deallocate(double ***, int, int, int);

  /// Update transformation matrix
  virtual void UpdateMatrix();

public:

  /// Returns the of control points in x
  virtual int GetX() const;

  /// Returns the of control points in y
  virtual int GetY() const;

  /// Returns the of control points in z
  virtual int GetZ() const;

  /// Returns the number of parameters of the transformation
  virtual int NumberOfDOFs() const;

  /// Get the control point spacing (in mm)
  virtual void GetSpacing(double &, double &, double &) const;

  /// Put orientation of free-form deformation
  virtual void  PutOrientation(double *, double *, double *);

  /// Get orientation of free-form deformation
  virtual void  GetOrientation(double *, double *, double *) const;

  /// Puts a control point value
  virtual void   Put(int, double);

  /// Gets a control point value
  virtual void   Put(int, int, int, double, double, double);

  /// Gets a control point value
  virtual double Get(int) const;

  /// Gets a control point value
  virtual void   Get(int, int, int, double &, double &, double &) const;

  /// Puts a control point status
  virtual void   PutStatus(int, int, int, _Status);

  /// Puts a control point status
  virtual void   PutStatus(int, int, int, _Status, _Status, _Status);

  /// Puts a control point status
  virtual void   PutStatus(int, _Status);

  /// Gets a control point status
  virtual void   GetStatus(int, int, int, _Status &, _Status &, _Status &);

  /// Gets a control point status
  virtual _Status GetStatus(int);

  /// Calculate the bending energy of the transformation
  virtual double Bending(double x, double y, double z, double t) = 0;

  /// Transforms world coordinates (in mm) to FFD coordinates
  virtual void WorldToLattice(double &, double &, double &) const;

  /// Transforms world coordinates (in mm) to FFD coordinates
  virtual void WorldToLattice(irtkPoint &) const;

  /// Transforms FFD coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(double &, double &, double &) const;

  /// Transforms FFD coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(irtkPoint &) const;

  /// Transforms index of control points to FFD coordinates
  virtual void IndexToLattice(int index, int& i, int& j, int& k) const;

  /// Transforms  FFD coordinates to index of control point
  virtual int  LatticeToIndex(int i, int j, int k) const;

  /// Returns the control point location (in mm)
  virtual void  ControlPointLocation(int, double &, double &, double &) const;

  /// Returns the control point location (in mm)
  virtual irtkPoint ControlPointLocation(int) const;

  /// Put the bounding box for FFD (in mm)
  virtual void PutBoundingBox(irtkPoint, irtkPoint);

  /// Put the bounding box for FFD (in mm)
  virtual void PutBoundingBox(double, double, double,
                              double, double, double);

  /// Returns the bounding box for FFD (in mm)
  virtual void BoundingBox(irtkPoint &, irtkPoint &) const;

  /// Returns the bounding box for FFD (in mm)
  virtual void BoundingBox(double &, double &, double &,
                           double &, double &, double &) const;

  /** Returns the bounding box for a control point (in mm). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(int, irtkPoint &, irtkPoint &, double = 1) const = 0;

  /** Returns the bounding box for a control point (in mm). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(int, double &, double &, double &,
                           double &, double &, double &, double = 1) const = 0;

  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(irtkGreyImage *, int, int &, int &, int &,
                           int &, int &, int &, double = 1) const = 0;

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *,
                             double *, double *, double *, int) = 0;

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dxs The x-displacements at each control point.
      \param dys The y-displacements at each control point.
      \param dzs The z-displacements at each control point. */
  virtual void Interpolate(double* dxs, double* dys, double* dzs) = 0;

  /// Inverts the transformation (abstract)
  virtual double Inverse(double &, double &, double &, double, double = 0.01);

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity();

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

};

inline int irtkFreeFormTransformation3D::GetX() const
{
  return _x;
}

inline int irtkFreeFormTransformation3D::GetY() const
{
  return _y;
}

inline int irtkFreeFormTransformation3D::GetZ() const
{
  return _z;
}

inline int irtkFreeFormTransformation3D::NumberOfDOFs() const
{
  return 3*_x*_y*_z;
}

inline double irtkFreeFormTransformation3D::Get(int index) const
{
  int i, j, k;

  if (index >= 3*_x*_y*_z) {
    cerr << "irtkFreeFormTransformation3D::Get: No such dof" << endl;
    exit(1);
  }
  if (index < _x*_y*_z) {
    i = index/(_y*_z);
    j = index%(_y*_z)/_z;
    k = index%(_y*_z)%_z;
    return _xdata[k][j][i];
  } else {
    if (index < 2*_x*_y*_z) {
      index -= _x*_y*_z;
      i = index/(_y*_z);
      j = index%(_y*_z)/_z;
      k = index%(_y*_z)%_z;
      return _ydata[k][j][i];
    } else {
      index -= 2*_x*_y*_z;
      i = index/(_y*_z);
      j = index%(_y*_z)/_z;
      k = index%(_y*_z)%_z;
      return _zdata[k][j][i];
    }
  }
}

inline void irtkFreeFormTransformation3D::Get(int i, int j, int k, double &x, double &y, double &z) const
{
  if ((i < 0) || (i >= _x) || (j < 0) || (j >= _y) || (k < 0) || (k >= _z)) {
    cerr << "irtkFreeFormTransformation3D::Get: No such dof" << endl;
    exit(1);
  }
  x = _xdata[k][j][i];
  y = _ydata[k][j][i];
  z = _zdata[k][j][i];
}

inline void irtkFreeFormTransformation3D::Put(int index, double x)
{
  irtkPoint p;
  int i, j, k;

  if (index >= 3*_x*_y*_z) {
    cerr << "irtkFreeFormTransformation3D::Put: No such dof" << endl;
    exit(1);
  }
  if (index < _x*_y*_z) {
    i = index/(_y*_z);
    j = index%(_y*_z)/_z;
    k = index%(_y*_z)%_z;
    _xdata[k][j][i] = x;
  } else {
    if (index < 2*_x*_y*_z) {
      index -= _x*_y*_z;
      i = index/(_y*_z);
      j = index%(_y*_z)/_z;
      k = index%(_y*_z)%_z;
      _ydata[k][j][i] = x;
    } else {
      index -= 2*_x*_y*_z;
      i = index/(_y*_z);
      j = index%(_y*_z)/_z;
      k = index%(_y*_z)%_z;
      _zdata[k][j][i] = x;
    }
  }
}

inline void irtkFreeFormTransformation3D::Put(int i, int j, int k, double x, double y, double z)
{
  if ((i < 0) || (i >= _x) || (j < 0) || (j >= _y) || (k < 0) || (k >= _z)) {
    cerr << "irtkFreeFormTransformation3D::Put: No such dof" << endl;
    exit(1);
  }
  _xdata[k][j][i] = x;
  _ydata[k][j][i] = y;
  _zdata[k][j][i] = z;
}

inline void irtkFreeFormTransformation3D::GetSpacing(double &dx, double &dy, double &dz) const
{
  dx = _dx;
  dy = _dy;
  dz = _dz;
}

inline void irtkFreeFormTransformation3D::PutOrientation(double *xaxis, double *yaxis, double *zaxis)
{
  _xaxis[0] = xaxis[0];
  _xaxis[1] = xaxis[1];
  _xaxis[2] = xaxis[2];
  _yaxis[0] = yaxis[0];
  _yaxis[1] = yaxis[1];
  _yaxis[2] = yaxis[2];
  _zaxis[0] = zaxis[0];
  _zaxis[1] = zaxis[1];
  _zaxis[2] = zaxis[2];

  this->UpdateMatrix();
}

inline void irtkFreeFormTransformation3D::GetOrientation(double *xaxis, double *yaxis, double *zaxis) const
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

inline void irtkFreeFormTransformation3D::PutBoundingBox(double x1, double y1, double z1, double x2, double y2, double z2)
{
  // Initialize control point domain
  _origin._x = (x2 + x1) / 2.0;
  _origin._y = (y2 + y1) / 2.0;
  _origin._z = (z2 + z1) / 2.0;

  double a = x1 * _xaxis[0] + y1 * _xaxis[1] + z1 * _xaxis[2];
  double b = x1 * _yaxis[0] + y1 * _yaxis[1] + z1 * _yaxis[2];
  double c = x1 * _zaxis[0] + y1 * _zaxis[1] + z1 * _zaxis[2];
  x1 = a;
  y1 = b;
  z1 = c;
  a = x2 * _xaxis[0] + y2 * _xaxis[1] + z2 * _xaxis[2];
  b = x2 * _yaxis[0] + y2 * _yaxis[1] + z2 * _yaxis[2];
  c = x2 * _zaxis[0] + y2 * _zaxis[1] + z2 * _zaxis[2];
  x2 = a;
  y2 = b;
  z2 = c;

  // Initialize control point spacing
  _dx = (x2 - x1) / (_x - 1);
  _dy = (y2 - y1) / (_y - 1);
  if (z2 > z1) {
    _dz = (z2 - z1) / (_z - 1);
  } else {
    _dz = 1;
  }

  this->UpdateMatrix();
}

inline void irtkFreeFormTransformation3D::PutBoundingBox(irtkPoint p1, irtkPoint p2)
{
  this->PutBoundingBox(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z);
}

inline void irtkFreeFormTransformation3D::WorldToLattice(double &x, double &y, double &z) const
{
  double a, b, c;

  // Pre-multiply point with transformation matrix
  a = _matW2L(0, 0)*x+_matW2L(0, 1)*y+_matW2L(0, 2)*z+_matW2L(0, 3);
  b = _matW2L(1, 0)*x+_matW2L(1, 1)*y+_matW2L(1, 2)*z+_matW2L(1, 3);
  c = _matW2L(2, 0)*x+_matW2L(2, 1)*y+_matW2L(2, 2)*z+_matW2L(2, 3);

  // Copy result back
  x = a;
  y = b;
  z = c;
}

inline void irtkFreeFormTransformation3D::WorldToLattice(irtkPoint &p) const
{
  this->WorldToLattice(p._x, p._y, p._z);
}

inline void irtkFreeFormTransformation3D::LatticeToWorld(double &x, double &y, double &z) const
{
  double a, b, c;

  // Pre-multiply point with transformation matrix
  a = _matL2W(0, 0)*x+_matL2W(0, 1)*y+_matL2W(0, 2)*z+_matL2W(0, 3);
  b = _matL2W(1, 0)*x+_matL2W(1, 1)*y+_matL2W(1, 2)*z+_matL2W(1, 3);
  c = _matL2W(2, 0)*x+_matL2W(2, 1)*y+_matL2W(2, 2)*z+_matL2W(2, 3);

  // Copy result back
  x = a;
  y = b;
  z = c;
}

inline void irtkFreeFormTransformation3D::LatticeToWorld(irtkPoint &p) const
{
  this->LatticeToWorld(p._x, p._y, p._z);
}

inline void irtkFreeFormTransformation3D::IndexToLattice(int index, int& i, int& j, int& k) const
{

  if (index >= _x*_y*_z) {
    index -= _x*_y*_z;
    if (index >= _x*_y*_z) {
      index -= _x*_y*_z;
    }
  }
  i = index/(_y*_z);
  j = index%(_y*_z)/_z;
  k = index%(_y*_z)%_z;
}

inline int irtkFreeFormTransformation3D::LatticeToIndex(int i, int j, int k) const
{
  return i * _y * _z + j * _z + k;
}

inline const char *irtkFreeFormTransformation3D::NameOfClass()
{
  return "irtkFreeFormTransformation3D";
}

#endif
