/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKPOINT_4D

#define _IRTKPOINT_4D

/**

  4D Point class.

*/

class irtkPoint4D : public irtkObject
{

public:

  /// x coordinate of Point
  double _x;

  /// y coordinate of Point
  double _y;

  /// z coordinate of Point
  double _z;

  /// t coordinate of Point
  double _t;

  //
  // Constructors and destructor
  //

  /// Constructor
  irtkPoint4D();

  /// Constructor with three coordinates
  irtkPoint4D(double, double, double,double);

  /// Constructor with Point
  irtkPoint4D(const irtkPoint4D &);

  /// Constructor with Vector
  irtkPoint4D(const irtkVector&);

  /// Default destructor
  virtual ~irtkPoint4D(void);

  //
  // Operators for Point
  //

  /// Copy operator for point
  irtkPoint4D& operator =(const irtkPoint4D&);

  /// Substraction operator for point
  irtkPoint4D& operator-=(const irtkPoint4D&);

  /// Addition operator for point
  irtkPoint4D& operator+=(const irtkPoint4D&);

  /// Multiplication operator for point
  irtkPoint4D& operator*=(const irtkPoint4D&);

  /// Division operator for point
  irtkPoint4D& operator/=(const irtkPoint4D&);

  /// Return result of point substraction
  irtkPoint4D  operator- (const irtkPoint4D&);

  /// Return result of point addition
  irtkPoint4D  operator+ (const irtkPoint4D&);

  /// Return result of point multiplication
  irtkPoint4D  operator* (const irtkPoint4D&);

  /// Return result of point division
  irtkPoint4D  operator/ (const irtkPoint4D&);

  //
  // Operators for comparison
  //

  /// Comparison operator ==
  int    operator==(const irtkPoint4D&);

  /// Comparison operator != (if USE_STL is defined, negate == operator)
  int    operator!=(const irtkPoint4D&);

  /// Comparison operator <
  int    operator<(const irtkPoint4D&);

  //
  // Operators for double
  //

  /// Substraction of double
  irtkPoint4D& operator-=(double);

  /// Addition of double
  irtkPoint4D& operator+=(double);

  /// Multiplication with double
  irtkPoint4D& operator*=(double);

  /// Division by double
  irtkPoint4D& operator/=(double);

  // Return result of substraction of double
  irtkPoint4D  operator- (double);

  // Return result of addition of double
  irtkPoint4D  operator+ (double);

  // Return result of multiplication with double
  irtkPoint4D  operator* (double);

  // Return result of division by double
  irtkPoint4D  operator/ (double);

  //
  // Operators for Vector
  //

  /// Copy operator for vectors
  irtkPoint4D& operator =(const irtkVector&);

  /// Substraction operator for vectors
  irtkPoint4D& operator-=(const irtkVector&);

  /// Addition operator for vectors
  irtkPoint4D& operator+=(const irtkVector&);

  /// Multiplication operator for vectors (componentwise)
  irtkPoint4D& operator*=(const irtkVector&);

  /// Division operator for vectors (componentwise)
  irtkPoint4D& operator/=(const irtkVector&);

  // Return result of vector substraction
  irtkPoint4D  operator- (const irtkVector&);

  // Return result of vector addition
  irtkPoint4D  operator+ (const irtkVector&);

  // Return result of vector multiplication
  irtkPoint4D  operator* (const irtkVector&);

  // Return result of vector division
  irtkPoint4D  operator/ (const irtkVector&);

  //
  // Operators for Matrix
  //

  /// Point multiplication operator for matrices
  irtkPoint4D& operator*=(const irtkMatrix&);

  /// Return result from Matrix multiplication
  irtkPoint4D  operator* (const irtkMatrix&);

  //
  // Distance methods
  //

  /// Distance from origin
  double  Distance(void) const;

  /// Distance from point
  double  Distance(const irtkPoint4D&) const;

  //
  // I/O methods
  //

  // Interface to output stream
  friend ostream& operator<< (ostream&, const irtkPoint4D&);

  /// Interface to input stream
  friend istream& operator>> (istream&, irtkPoint4D&);
};

inline irtkPoint4D::irtkPoint4D()
{
  _x = 0;
  _y = 0;
  _z = 0;
  _t = 0;
}

inline irtkPoint4D::irtkPoint4D(double x, double y, double z,double t)
{
  _x = x;
  _y = y;
  _z = z;
  _t = t;
}

inline irtkPoint4D::irtkPoint4D(const irtkPoint4D& p)
{
  _x = p._x;
  _y = p._y;
  _z = p._z;
  _t = p._t;
}

inline irtkPoint4D::irtkPoint4D(const irtkVector& v)
{
  if ((v.Rows() < 0) || (v.Rows() > 4)) {
    cerr << "irtkPoint4D::irtkPoint4D(const irtkVector&) Illegal dimension: " << v.Rows() << endl;
    exit(1);
  } else {
    if (v.Rows() == 1) {
      _x = v(0);
      _y = 0;
      _z = 0;
      _t = 0;
    }
    if (v.Rows() == 2) {
      _x = v(0);
      _y = v(1);
      _z = 0;
      _t = 0;
    }
    if (v.Rows() == 3) {
      _x = v(0);
      _y = v(1);
      _z = v(2);
      _t = 0;
    }
    if (v.Rows() == 4) {
      _x = v(0);
      _y = v(1);
      _z = v(2);
      _t = v(4);
    }
  }
}

inline irtkPoint4D::~irtkPoint4D()
{
}

inline irtkPoint4D& irtkPoint4D::operator =(const irtkPoint4D& p)
{
  _x = p._x;
  _y = p._y;
  _z = p._z;
  _t = p._t;
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator+=(const irtkPoint4D& p)
{
  _x += p._x;
  _y += p._y;
  _z += p._z;
  _t += p._t;
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator-=(const irtkPoint4D& p)
{
  _x -= p._x;
  _y -= p._y;
  _z -= p._z;
  _t -= p._t;
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator*=(const irtkPoint4D& p)
{
  _x *= p._x;
  _y *= p._y;
  _z *= p._z;
  _t *= p._t;
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator/=(const irtkPoint4D& p)
{
  _x /= p._x;
  _y /= p._y;
  _z /= p._z;
  _t /= p._t;
  return *this;
}

inline irtkPoint4D irtkPoint4D::operator+ (const irtkPoint4D& p)
{
  irtkPoint4D tmp;

  tmp._x = _x + p._x;
  tmp._y = _y + p._y;
  tmp._z = _z + p._z;
  tmp._t = _t + p._t;
  return tmp;
}

inline irtkPoint4D irtkPoint4D::operator- (const irtkPoint4D& p)
{
  irtkPoint4D tmp;

  tmp._x = _x - p._x;
  tmp._y = _y - p._y;
  tmp._z = _z - p._z;
  tmp._t = _t - p._t;
  return tmp;
}

inline irtkPoint4D irtkPoint4D::operator* (const irtkPoint4D& p)
{
  irtkPoint4D tmp;

  tmp._x = _x * p._x;
  tmp._y = _y * p._y;
  tmp._z = _z * p._z;
  tmp._t = _t * p._t;
  return tmp;
}

inline irtkPoint4D irtkPoint4D::operator/ (const irtkPoint4D& p)
{
  irtkPoint4D tmp;

  tmp._x = _x / p._x;
  tmp._y = _y / p._y;
  tmp._z = _z / p._z;
  tmp._t = _t / p._t;
  return tmp;
}

inline int irtkPoint4D::operator==(irtkPoint4D const &p)
{
  return ((_x == p._x) && (_y == p._y) && (_z == p._z) && (_t==p._t) );
}

inline int irtkPoint4D::operator!=(irtkPoint4D const &p)
{
  return ((_x != p._x) || (_y != p._y) || (_z != p._z) || (_t!=p._t) );
}

inline irtkPoint4D& irtkPoint4D::operator+=(double x)
{
  _x += x;
  _y += x;
  _z += x;
  _t += x;
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator-=(double x)
{
  _x -= x;
  _y -= x;
  _z -= x;
  _t -= x;
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator/=(double x)
{
  _x /= x;
  _y /= x;
  _z /= x;
  _t /= x;
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator*=(double x)
{
  _x *= x;
  _y *= x;
  _z *= x;
  _t *= x;
  return *this;
}

inline irtkPoint4D irtkPoint4D::operator+ (double x)
{
  irtkPoint4D p;

  p._x = _x + x;
  p._y = _y + x;
  p._z = _z + x;
  p._t = _t + x;
  return p;
}

inline irtkPoint4D irtkPoint4D::operator- (double x)
{
  irtkPoint4D p;

  p._x = _x - x;
  p._y = _y - x;
  p._z = _z - x;
  p._t = _t - x;
  return p;
}

inline irtkPoint4D irtkPoint4D::operator* (double x)
{
  irtkPoint4D p;

  p._x = _x * x;
  p._y = _y * x;
  p._z = _z * x;
  p._t = _t * x;
  return p;
}

inline irtkPoint4D irtkPoint4D::operator/ (double x)
{
  irtkPoint4D p;

  p._x = _x / x;
  p._y = _y / x;
  p._z = _z / x;
  p._t = _t / x;
  return p;
}

inline irtkPoint4D& irtkPoint4D::operator+=(const irtkVector& v)
{
  if (v.Rows() != 4) {
    cerr << "irtkPoint4D::operator+=(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  _x += v(0);
  _y += v(1);
  _z += v(2);
  _t += v(3);
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator-=(const irtkVector& v)
{
  if (v.Rows() != 4) {
    cerr << "irtkPoint4D::operator-=(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  _x -= v(0);
  _y -= v(1);
  _z -= v(2);
  _t -= v(3);
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator*=(const irtkVector& v)
{
  if (v.Rows() != 4) {
    cerr << "irtkPoint4D::operator*=(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  _x *= v(0);
  _y *= v(1);
  _z *= v(2);
  _t *= v(3);
  return *this;
}

inline irtkPoint4D& irtkPoint4D::operator/=(const irtkVector& v)
{
  if (v.Rows() != 4) {
    cerr << "irtkPoint4D::operator/=(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  _x /= v(0);
  _y /= v(1);
  _z /= v(2);
  _t /= v(3);
  return *this;
}

inline irtkPoint4D irtkPoint4D::operator+ (const irtkVector& v)
{
  irtkPoint4D tmp;

  if (v.Rows() != 4) {
    cerr << "irtkPoint4D::operator+(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x + v(0);
  tmp._y = _y + v(1);
  tmp._z = _z + v(2);
  tmp._t = _t + v(3);
  return tmp;
}

inline irtkPoint4D irtkPoint4D::operator- (const irtkVector& v)
{
  irtkPoint4D tmp;

  if (v.Rows() != 4) {
    cerr << "irtkPoint4D::operator-(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x - v(0);
  tmp._y = _y - v(1);
  tmp._z = _z - v(2);
  tmp._t = _t - v(3);
  return tmp;
}

inline irtkPoint4D irtkPoint4D::operator* (const irtkVector& v)
{
  irtkPoint4D tmp;

  if (v.Rows() != 4) {
    cerr << "irtkPoint4D::operator*(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x * v(0);
  tmp._y = _y * v(1);
  tmp._z = _z * v(2);
  tmp._t = _t * v(3);
  return tmp;
}

inline irtkPoint4D irtkPoint4D::operator/ (const irtkVector& v)
{
  irtkPoint4D tmp;

  if (v.Rows() != 4) {
    cerr << "irtkPoint4D::operator/(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x / v(0);
  tmp._y = _y / v(1);
  tmp._z = _z / v(2);
  tmp._t = _t / v(3);
  return tmp;
}

inline irtkPoint4D irtkPoint4D::operator* (const irtkMatrix& m)
{
  irtkPoint4D tmp;

  if ((m.Rows() != 5) && (m.Cols() != 5)) {
    cerr << "irtkPoint4D::operator*(const irtkMatrix& m): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x * m(0, 0) + _y * m(0, 1) + _z * m(0, 2) + _t * m(0, 3) + m(0,4);
  tmp._y = _x * m(1, 0) + _y * m(1, 1) + _z * m(1, 2) + _t * m(1, 3) + m(1,4);
  tmp._z = _x * m(2, 0) + _y * m(2, 1) + _z * m(2, 2) + _t * m(2, 3) + m(2,4);
  tmp._t = _x * m(3, 0) + _y * m(3, 1) + _z * m(3, 2) + _t * m(3, 3) + m(3,4);

  return tmp;
}

inline irtkPoint4D& irtkPoint4D::operator*=(const irtkMatrix& m)
{
  irtkPoint4D tmp;

  if ((m.Rows() != 4) && (m.Cols() != 4)) {
    cerr << "irtkPoint4D::operator*(const irtkMatrix& m): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x * m(0, 0) + _y * m(0, 1) + _z * m(0, 2) + _t * m(0, 3) + m(0,4);
  tmp._y = _x * m(1, 0) + _y * m(1, 1) + _z * m(1, 2) + _t * m(1, 3) + m(1,4);
  tmp._z = _x * m(2, 0) + _y * m(2, 1) + _z * m(2, 2) + _t * m(2, 3) + m(2,4);
  tmp._t = _x * m(3, 0) + _y * m(3, 1) + _z * m(3, 2) + _t * m(3, 3) + m(3,4);
  *this  = tmp;

  return *this;
}

inline double irtkPoint4D::Distance(void) const
{
  return sqrt(_x*_x + _y*_y  + _z *_z + _t*_t);
}

inline double irtkPoint4D::Distance(const irtkPoint4D &p) const
{
  return sqrt((_x - p._x)*(_x - p._x) + (_y - p._y)*(_y - p._y) +
              (_z - p._z)*(_z - p._z) + (_t - p._t)*(_t - p._t));
}

#endif
