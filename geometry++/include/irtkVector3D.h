/*=========================================================================
 
  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$
 
=========================================================================*/

#ifndef _IRTKVECTOR3D_H

#define _IRTKVECTOR3D_H

#include <cmath>
#include <limits>

#include <irtkObject.h>

/** Represents a 3D vector. */
template <typename T> class irtkVector3D : public irtkObject
{
public:
  /** The x-component. */
  T _x;

  /** The y-component. */
  T _y;

  /** The z-component. */
  T _z;

  /** Constructor. */
  irtkVector3D();

  /** Constructor. */
  irtkVector3D(T x, T y, T z);

  /** Copy Constructor. */
  irtkVector3D(const irtkVector3D &);

  /** Operator for dividing one vector by another. */
  irtkVector3D& operator=(const irtkVector3D& v);

  /** Operator for multiplying by a scalar. */
  irtkVector3D operator*(double s);

  /** Operator for adding two vectors. */
  irtkVector3D operator+(const irtkVector3D& v);

  /** Operator for subtraction. */
  irtkVector3D operator-(const irtkVector3D& v);

  /** Operator for multiplying by a scalar. */
  irtkVector3D& operator*=(double s);

  /** Operator for mulityplying by a vector. */
  irtkVector3D& operator*=(const irtkVector3D& v);

  /** Operator for adding a vector. */
  irtkVector3D& operator+=(const irtkVector3D& v);

  /** Operator for subtracting a vector. */
  irtkVector3D& operator-=(const irtkVector3D& v);

  /** Operator for testing equality of two vectors. */
  bool operator==(const irtkVector3D& v);

  /** Operator for testing non-equality of two vector. */
  bool operator!=(const irtkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator<(const irtkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator>(const irtkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator<=(const irtkVector3D& v);

  /** Operator for comparing sizes of vectors. */
  bool operator>=(const irtkVector3D& v);

  /** Operator for dividing one vector by another. */
  irtkVector3D& operator/=(const irtkVector3D& v);

  /** Operator for dividing one vector by another. */
  irtkVector3D operator/(const irtkVector3D& v);

  /** Operator for dividing a vector by a scalar. */
  irtkVector3D& operator/=(double s);

  /** Operator for dividing a vector by a scalar. */
  irtkVector3D operator/(double s);

  /** Normalizes the vector. */
  void Normalize();

  /** Returns the length of the vector. */
  double Length() const;

  /** Takes the cross-product of two vectors. */
  static irtkVector3D CrossProduct(const irtkVector3D& v1, const irtkVector3D& v2);

  /** Takes the dot-product of two vectors. */
  static double DotProduct(const irtkVector3D& v1, const irtkVector3D& v2);

};

template <typename T> inline irtkVector3D<T>::irtkVector3D()
{
  _x = 0;
  _y = 0;
  _z = 0;
}

template <typename T> inline irtkVector3D<T>::irtkVector3D(T x, T y, T z)
{
  _x = x;
  _y = y;
  _z = z;
}
template <typename T> inline irtkVector3D<T>::irtkVector3D(const irtkVector3D &v)
{
  _x = v._x;
  _y = v._y;
  _z = v._z;
}

template <typename T> inline irtkVector3D<T>& irtkVector3D<T>::operator=(const irtkVector3D<T>& v)
{
  _x = v._x;
  _y = v._y;
  _z = v._z;
  return *this;
}

template <typename T> inline irtkVector3D<T> irtkVector3D<T>::operator*(double s)
{
  irtkVector3D<T> r;

  r._x = static_cast<T>(_x*s);
  r._y = static_cast<T>(_y*s);
  r._z = static_cast<T>(_z*s);

  return r;
}

template <typename T> inline irtkVector3D<T> irtkVector3D<T>::operator+(const irtkVector3D<T>& v)
{
  irtkVector3D<T> r;

  r._x = _x + v._x;
  r._y = _y + v._y;
  r._z = _z + v._z;

  return r;
}

template <typename T> inline irtkVector3D<T> irtkVector3D<T>::operator-(const irtkVector3D<T>& v)
{
  irtkVector3D<T> r;

  r._x = _x - v._x;
  r._y = _y - v._y;
  r._z = _z - v._z;

  return r;
}


template <typename T> inline irtkVector3D<T>& irtkVector3D<T>::operator*=(double s)
{
  _x = static_cast<T>(_x*s);
  _y = static_cast<T>(_y*s);
  _z = static_cast<T>(_z*s);

  return *this;
}

template <typename T> inline irtkVector3D<T>& irtkVector3D<T>::operator*=(const irtkVector3D<T>& v)
{
  _x *= v._x;
  _y *= v._y;
  _z *= v._z;

  return *this;
}

template <typename T> inline irtkVector3D<T>& irtkVector3D<T>::operator+=(const irtkVector3D<T>& v)
{
  _x += v._x;
  _y += v._y;
  _z += v._z;

  return *this;
}

template <typename T> inline irtkVector3D<T>& irtkVector3D<T>::operator-=(const irtkVector3D<T>& v)
{
  _x -= v._x;
  _y -= v._y;
  _z -= v._z;

  return *this;
}

template <typename T> inline bool irtkVector3D<T>::operator==(const irtkVector3D<T>& v)
{
  return ((_z == v._z) && (_y == v._y) && (_x == v._x));
}

template <typename T> inline bool irtkVector3D<T>::operator!=(const irtkVector3D<T>& v)
{
  return ((_z != v._z) || (_y != v._y) || (_x != v._x));
}

template <typename T> inline bool irtkVector3D<T>::operator<(const irtkVector3D<T>& v)
{
  return ((_z < v._z) ||
          ((_z == v._z) && (_y < v._y)) ||
          ((_z == v._z) && (_y == v._y) && (_x < v._x)));
}

template <typename T> inline bool irtkVector3D<T>::operator>(const irtkVector3D<T>& v)
{
  return ((_z > v._z) ||
          ((_z == v._z) && (_y > v._y)) ||
          ((_z == v._z) && (_y == v._y) && (_x > v._x)));
}

template <typename T> inline bool irtkVector3D<T>::operator<=(const irtkVector3D<T>& v)
{
  return ((*this < v) || (*this == v));
}

template <typename T> inline bool irtkVector3D<T>::operator>=(const irtkVector3D<T>& v)
{
  return ((*this > v) || (*this == v));
}

template <typename T> inline irtkVector3D<T>& irtkVector3D<T>::operator/=(double s)
{
  _x = static_cast<T>(_x/s);
  _y = static_cast<T>(_y/s);
  _z = static_cast<T>(_z/s);

  return *this;
}

template <typename T> inline irtkVector3D<T> irtkVector3D<T>::operator/(double s)
{
  return irtkVector3D<T>(static_cast<T>(_x/s), static_cast<T>(_y/s), static_cast<T>(_z/s));
}

template <typename T> inline void irtkVector3D<T>::Normalize()
{
  double length = sqrt(static_cast<double>(_x*_x + _y*_y + _z*_z));

  if (length != 0) {
    _x = static_cast<T>(_x/length);
    _y = static_cast<T>(_y/length);
    _z = static_cast<T>(_z/length);
  }
}

template <typename T> inline double irtkVector3D<T>::Length() const
{
  return sqrt(static_cast<double>(_x*_x + _y*_y + _z*_z));
}

template<typename T> inline irtkVector3D<T> irtkVector3D<T>::CrossProduct(const irtkVector3D<T>& v1, const irtkVector3D<T>& v2)
{
  irtkVector3D<T> cp;

  cp._x = v1._y*v2._z - v1._z*v2._y;
  cp._y = v1._z*v2._x - v1._x*v2._z;
  cp._z = v1._x*v2._y - v1._y*v2._x;

  return cp;
}

template <typename T> inline double irtkVector3D<T>::DotProduct(const irtkVector3D<T>& v1, const irtkVector3D<T>& v2)
{
  return v1._x*v2._x + v1._y*v2._y + v1._z*v2._z;
}

#endif // __IRTKVECTOR3D_H
