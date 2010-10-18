/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKQUATERNION_H

#define _IRTKQUATERNION_H

#include <irtkCommon.h>
#include <irtkGeometry.h>

/** Represents a quaternion. */
class irtkQuaternion : public irtkObject
{
public:
  /** The real part of the quaternion. */
  double _t;

  /** The first imaginary part. */
  double _u;

  /** The second imaginary part. */
  double _v;

  /** The third imaginary part. */
  double _w;

  /** Default constructor creates a quaternion with components (0, 0, 0, 0). */
  irtkQuaternion();

  /** Constructor creates a quaternion with components (t, u, v, w). */
  irtkQuaternion(double t, double u, double v, double w);

  /** Copy constructor. */
  irtkQuaternion(const irtkQuaternion& q);

  /** Assignment operator. */
  irtkQuaternion& operator=(const irtkQuaternion& q);

  /** Adds two quaternions. */
  friend irtkQuaternion operator+(const irtkQuaternion& a, const irtkQuaternion& b);

  /** Subtracts two quaternions. */
  friend irtkQuaternion operator-(const irtkQuaternion& a, const irtkQuaternion& b);

  /** Multiplies two quaternions. */
  friend irtkQuaternion operator*(const irtkQuaternion& a, const irtkQuaternion& b);

  /** Multiplies a quaternion by a real number. */
  friend irtkQuaternion operator*(const irtkQuaternion& a, double b);

  /** Outputs the quaternion. */
  friend ostream& operator<<(ostream& os, const irtkQuaternion& q);

  /** Returns the conjugate of a quaternion. */
  irtkQuaternion Conjugate() const;

  /** Returns the inverse of the quaternion. */
  irtkQuaternion Inverse() const;

  /** Normalizes the quaternion. */
  void Normalize();

  /** Returns the length of the quaternion. */
  double Length() const;

  /** Computes the rotation matrix corresponding to this quaternion. Note
      that the matrix is a 4x4 homogeneous matrix. */
  irtkMatrix ComputeRotationMatrix() const;

  /** Computes the rotation matrix needed to a rotate a point by an angle
      alpha about the axis with direction (u, v, w) passing through the
      point (ox, oy, oz). The direction of rotation is the direction a
      right-handed corkscrew needs to be rotated in order to move along the
      direction of the axis. Note that the matrix returned is a 4x4
      homogeneous matrix. */
  static irtkMatrix ComputeTransformationMatrix(double alpha,
      double u, double v, double w, double ox, double oy, double oz);

  /** Computes the rotation matrix needed to rotate a point by an angle alpha
      about the axis with direction (u, v, w) passing through the origin.
      Direction of rotation is direction a right-handed corkscrew needs to be
      rotated in order to move along direction of axis. Note that the
      rotation matrix is a 4x4 homogeneous matrix. */
  static irtkMatrix ComputeRotationMatrix(double alpha, double u, double v, double w);

  /** Computes the quaternion which represents a rotation by an angle alpha
      about the axis with direction (u, v, w) passing through the origin.
      Direction of rotation is direction a right-handed corkscrew needs to be
      rotated in order to move along direction of axis. */
  static irtkQuaternion ComputeRotationQuaternion(double alpha, double u, double v, double w);

  /** Returns the name of the class. */
  virtual const char* NameOfClass();
};

inline irtkQuaternion::irtkQuaternion()
{
  _t = 0;
  _u = 0;
  _v = 0;
  _w = 0;
}

inline irtkQuaternion::irtkQuaternion(double t, double u, double v, double w)
{
  _t = t;
  _u = u;
  _v = v;
  _w = w;
}

inline irtkQuaternion::irtkQuaternion(const irtkQuaternion& q) : irtkObject(q)
{
  _t = q._t;
  _u = q._u;
  _v = q._v;
  _w = q._w;
}

inline irtkQuaternion& irtkQuaternion::operator=(const irtkQuaternion& q)
{
  if (&q != this) {
    _t = q._t;
    _u = q._u;
    _v = q._v;
    _w = q._w;
  }

  return *this;
}

inline irtkQuaternion operator+(const irtkQuaternion& a, const irtkQuaternion& b)
{
  irtkQuaternion q;

  q._t = a._t + b._t;
  q._u = a._u + b._u;
  q._v = a._v + b._v;
  q._w = a._w + b._w;

  return q;
}

inline irtkQuaternion operator-(const irtkQuaternion& a, const irtkQuaternion& b)
{
  irtkQuaternion q;

  q._t = a._t - b._t;
  q._u = a._u - b._u;
  q._v = a._v - b._v;
  q._w = a._w - b._w;

  return q;
}

inline irtkQuaternion operator*(const irtkQuaternion& a, const irtkQuaternion& b)
{
  irtkQuaternion q;

  q._t = a._t*b._t - a._u*b._u - a._v*b._v - a._w*b._w;
  q._u = a._t*b._u + a._u*b._t + a._v*b._w - a._w*b._v;
  q._v = a._t*b._v - a._u*b._w + a._v*b._t + a._w*b._u;
  q._w = a._t*b._w + a._u*b._v - a._v*b._u + a._w*b._t;

  return q;
}

inline irtkQuaternion operator*(const irtkQuaternion& a, double b)
{
  irtkQuaternion q;

  q._t = a._t*b;
  q._u = a._u*b;
  q._v = a._v*b;
  q._w = a._w*b;

  return q;
}

inline ostream& operator<<(ostream& os, const irtkQuaternion& q)
{
  os << "_t = " << q._t << " _u = " << q._u << " _v = " << q._v << " _w = " << q._w;

  return os;
}

inline irtkQuaternion irtkQuaternion::Conjugate() const
{
  irtkQuaternion q;

  q._t = _t;
  q._u = -_u;
  q._v = -_v;
  q._w = -_w;

  return q;
}

inline irtkQuaternion irtkQuaternion::Inverse() const
{
  irtkQuaternion q = Conjugate();

  double lenq = q.Length();
  q = q*(1.0/(lenq*lenq));

  return q;
}

inline void irtkQuaternion::Normalize()
{
  double length = Length();

  if (length > 0) {
    _t /= length;
    _u /= length;
    _v /= length;
    _w /= length;
  }
}

inline double irtkQuaternion::Length() const
{
  return sqrt(_t*_t + _u*_u + _v*_v + _w*_w);
}

inline const char* irtkQuaternion::NameOfClass()
{
  return "irtkQuaternion";
}

#endif // __IRTKQUATERNION_H
