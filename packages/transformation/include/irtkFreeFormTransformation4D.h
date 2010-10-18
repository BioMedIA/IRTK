/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKFREEFORMTRANSFORMATION4D_H

#define _IRTKFREEFORMTRANSFORMATION4D_H

#include <irtkGeometry.h>

#define FFDLOOKUPTABLESIZE 1000

/**
 * Class for 4D free form transformations
 *
 * This class implements a 4D free form transformation
 *
 */

class irtkFreeFormTransformation4D : public irtkFreeFormTransformation
{

protected:

  /// Number of control points in x
  int _x;

  /// Number of control points in y
  int _y;

  /// Number of control points in z
  int _z;

  /// Number of control points in t
  int _t;

  /// Spacing of control points in x (in mm)
  double _dx;

  /// Spacing of control points in y (in mm)
  double _dy;

  /// Spacing of control points in z (in mm)
  double _dz;

  /// Spacing of control points in t (in ms)
  double _dt;

  /// The minimum time value for which the transformation is defined.
  double _tMin;

  /// The maximum time value for which the transformation is defined.
  double _tMax;

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
  double ****_xdata;

  /// Displacement in the y-direction at the control points (in mm)
  double ****_ydata;

  /// Displacement in the z-direction at the control points (in mm)
  double ****_zdata;

  /// Allocate memory for control points
  static double ****Allocate  (double ****, int, int, int, int);

  /// Deallocate memory for control points
  static double ****Deallocate(double ****, int, int, int, int);

  /// Update transformation matrix
  virtual void UpdateMatrix();

public:

  /// Subdivide FFD
  virtual void Subdivide() = 0;

  /// Returns the of control points in x
  virtual int GetX() const;

  /// Returns the of control points in y
  virtual int GetY() const;

  /// Returns the of control points in z
  virtual int GetZ() const;

  /// Returns the of control points in t
  virtual int GetT() const;

  /// Returns the number of parameters of the transformation
  virtual int NumberOfDOFs() const;

  /// Get the control point spacing (in mm and ms)
  virtual void GetSpacing(double &, double &, double &, double &) const;

  /// Put orientation of free-form deformation
  virtual void  PutOrientation(double *, double *, double *);

  /// Get orientation of free-form deformation
  virtual void  GetOrientation(double *, double *, double *) const;

  /// Puts a control point value
  virtual void   Put(int, double);

  /// Gets a control point value
  virtual void   Put(int, int, int, int, double, double, double);

  /// Gets a control point value
  virtual double Get(int) const;

  /// Gets a control point value
  virtual void   Get(int, int, int, int, double &, double &, double &) const;

  /// Calculate the bending energy of the transformation
  virtual double Bending(double, double, double, double) = 0;

  /// Puts a control point status
  virtual void   PutStatus(int, int, int, int, _Status);

  /// Puts a control point status
  virtual void   PutStatus(int, int, int, int, _Status, _Status, _Status);

  /// Puts a control point status
  virtual void   PutStatus(int, _Status);

  /// Gets a control point status
  virtual void   GetStatus(int, int, int, int, _Status &, _Status &, _Status &);

  /// Gets a control point status
  virtual _Status GetStatus(int);

  /// Transforms world coordinates (in mm) to FFD coordinates
  virtual void WorldToLattice(double &, double &, double &) const;

  /// Transforms world coordinates (in mm) to FFD coordinates
  virtual void WorldToLattice(irtkPoint &) const;

  /// Transforms time coordinate (in ms) to FFD coordinates
  virtual double TimeToLattice(double) const;

  /// Transforms FFD coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(double &, double &, double &) const;

  /// Transforms FFD coordinates to world coordinates (in mm)
  virtual void LatticeToWorld(irtkPoint &) const;

  /// Transforms FFD coordinates to time coordinate (in ms)
  virtual double LatticeToTime(double) const;

  /// Transforms index of control points to FFD coordinates
  virtual void IndexToLattice(int index, int& i, int& j, int& k, int& l) const;

  /// Transforms  FFD coordinates to index of control point
  virtual int  LatticeToIndex(int i, int j, int k, int l) const;

  /// Returns the control point location (in mm)
  virtual void  ControlPointLocation(int, double &, double &, double &) const;

  /// Returns the control point location (in mm)
  virtual irtkPoint ControlPointLocation(int) const;

  /// Returns the bounding box for FFD (in mm)
  virtual void BoundingBox(irtkPoint &, irtkPoint &) const;

  /// Returns the bounding box for FFD (in mm)
  virtual void BoundingBox(double &, double &, double &,
                           double &, double &, double &) const;

  /** Returns the bounding box for a control point. The last parameter
   *  specifies what fraction of the bounding box to return. The default
   *  is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(int, double &, double &, double &, double &,
                           double &, double &, double &, double &, double = 1) const = 0;

  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(irtkGreyImage *, int, int &, int &, int &, int &,
                           int &, int &, int &, int &, double = 1) const = 0;

  /// Inverts the transformation (abstract)
  virtual double Inverse(double &, double &, double &, double = 0, double = 0.01);

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity();

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

};

inline int irtkFreeFormTransformation4D::GetX() const
{
  return _x;
}

inline int irtkFreeFormTransformation4D::GetY() const
{
  return _y;
}

inline int irtkFreeFormTransformation4D::GetZ() const
{
  return _z;
}

inline int irtkFreeFormTransformation4D::GetT() const
{
  return _t;
}

inline int irtkFreeFormTransformation4D::NumberOfDOFs() const
{
  return 3*_x*_y*_z*_t;
}

inline double irtkFreeFormTransformation4D::Get(int index) const
{
  int i, j, k, l;

  if (index >= 3*_x*_y*_z*_t) {
    cerr << "irtkFreeFormTransformation4D::Get: No such dof" << endl;
    exit(1);
  }

  if (index < _x*_y*_z*_t) {
    l = index/(_z*_y*_x);
    k = index%(_z*_y*_x)/(_y*_x);
    j = index%(_z*_y*_x)%(_y*_x)/_x;
    i = index%(_z*_y*_x)%(_y*_x)%_x;
    cout << "x = " << i << " " << j << " " << k << " " << l << endl;

    return _xdata[l][k][j][i];
  } else if (index < 2*_x*_y*_z*_t) {
    index -= _x*_y*_z*_t;

    l = index/(_z*_y*_x);
    k = index%(_z*_y*_x)/(_y*_x);
    j = index%(_z*_y*_x)%(_y*_x)/_x;
    i = index%(_z*_y*_x)%(_y*_x)%_x;
    cout << "y  = " << i << " " << j << " " << k << " " << l << endl;

    return _ydata[l][k][j][i];
  } else {
    index -= 2*_x*_y*_z*_t;

    l = index/(_z*_y*_x);
    k = index%(_z*_y*_x)/(_y*_x);
    j = index%(_z*_y*_x)%(_y*_x)/_x;
    i = index%(_z*_y*_x)%(_y*_x)%_x;
    cout << "z = " << i << " " << j << " " << k << " " << l << endl;

    return _zdata[l][k][j][i];
  }
}

inline void irtkFreeFormTransformation4D::Get(int i, int j, int k, int l, double &x, double &y, double &z) const
{
  if ((i < 0) || (i >= _x) || (j < 0) || (j >= _y) || (k < 0) || (k >= _z) || (l < 0) || (l >= _t)) {
    cerr << "irtkFreeFormTransformation4D::Get: No such dof" << endl;
    exit(1);
  }
  x = _xdata[l][k][j][i];
  y = _ydata[l][k][j][i];
  z = _zdata[l][k][j][i];
}

inline void irtkFreeFormTransformation4D::Put(int index, double x)
{
  int i, j, k, l;

  if (index >= 3*_x*_y*_z*_t) {
    cerr << "irtkFreeFormTransformation4D::Put: No such dof" << endl;
    exit(1);
  }

  if (index < _x*_y*_z*_t) {
    l = index/(_z*_y*_x);
    k = index%(_z*_y*_x)/(_y*_x);
    j = index%(_z*_y*_x)%(_y*_x)/_x;
    i = index%(_z*_y*_x)%(_y*_x)%_x;

    _xdata[l][k][j][i] = x;
  } else if (index < 2*_x*_y*_z*_t) {
    index -= _x*_y*_z*_t;

    l = index/(_z*_y*_x);
    k = index%(_z*_y*_x)/(_y*_x);
    j = index%(_z*_y*_x)%(_y*_x)/_x;
    i = index%(_z*_y*_x)%(_y*_x)%_x;

    _ydata[l][k][j][i] = x;
  } else {
    index -= 2*_x*_y*_z*_t;

    l = index/(_z*_y*_x) ;
    k = index%(_z*_y*_x)/(_y*_x) ;
    j = index%(_z*_y*_x)%(_y*_x)/_x;
    i = index%(_z*_y*_x)%(_y*_x)%_x;

    _zdata[l][k][j][i] = x;
  }
}

inline void irtkFreeFormTransformation4D::Put(int i, int j, int k, int l, double x, double y, double z)
{
  if ((i < 0) || (i >= _x) || (j < 0) || (j >= _y) || (k < 0) || (k >= _z) || (l < 0) || (l >= _t)) {
    cerr << "irtkFreeFormTransformation4D::Put: No such dof" << endl;
    exit(1);
  }
  _xdata[l][k][j][i] = x;
  _ydata[l][k][j][i] = y;
  _zdata[l][k][j][i] = z;
}

inline void irtkFreeFormTransformation4D::GetSpacing(double &dx, double &dy, double &dz, double &dt) const
{
  dx = _dx;
  dy = _dy;
  dz = _dz;
  dt = _dt;
}

inline void irtkFreeFormTransformation4D::PutOrientation(double *xaxis, double *yaxis, double *zaxis)
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

inline void irtkFreeFormTransformation4D::GetOrientation(double *xaxis, double *yaxis, double *zaxis) const
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

inline void irtkFreeFormTransformation4D::WorldToLattice(double &x, double &y, double &z) const
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

inline void irtkFreeFormTransformation4D::WorldToLattice(irtkPoint &p) const
{
  this->WorldToLattice(p._x, p._y, p._z);
}

inline double irtkFreeFormTransformation4D::TimeToLattice(double t) const
{
  return (t - _tMin)*(_t - 1)/(_tMax - _tMin);
}

inline void irtkFreeFormTransformation4D::LatticeToWorld(double &x, double &y, double &z) const
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

inline void irtkFreeFormTransformation4D::LatticeToWorld(irtkPoint &p) const
{
  this->LatticeToWorld(p._x, p._y, p._z);
}

inline double irtkFreeFormTransformation4D::LatticeToTime(double t) const
{
  return t*(_tMax - _tMin)/double(_t - 1)+_tMin;
}

inline void irtkFreeFormTransformation4D::IndexToLattice(int index, int& i, int& j, int& k, int& l) const
{
  if (index > _x*_y*_z*_t) {
    index -= _x*_y*_z*_t;
    if (index > _x*_y*_z*_t)
      index -= _x*_y*_z*_t;
  }

  l = index/(_z*_y*_x);
  k = index%(_z*_y*_x)/(_y*_x);
  j = index%(_z*_y*_x)%(_y*_x)/_x;
  i = index%(_z*_y*_x)%(_y*_x)%_x;
}

inline int irtkFreeFormTransformation4D::LatticeToIndex(int i, int j, int k, int l) const
{
  return l*_z*_y*_x + k*_y*_x + j*_x + i;
}

inline const char *irtkFreeFormTransformation4D::NameOfClass()
{
  return "irtkFreeFormTransformation4D";
}

#endif
