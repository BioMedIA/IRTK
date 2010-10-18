/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKVECTOR_H

#define _IRTKVECTOR_H

#ifdef USE_VXL
// VXL header file was loaded in irtkMatrix.h
#include <vnl/vnl_diag_matrix.h>
#endif

/**

  Vector class.

*/

class irtkVector : public irtkObject
{

protected:

  /// Number of rows
  int _rows;

  /// Data
  double *_vector;

public:

  /// Default constructor
  irtkVector();

  /// Constructor for given row dimensions
  irtkVector(int);

  /// Copy constructor
  irtkVector(const irtkVector &);

  /// Destructor
  ~irtkVector();

  /// Intialize matrix with number of rows
  void Initialize(int);

  //
  // Vector access functions
  //

  /// Returns number of rows
  int Rows() const;

  /// Puts vector value
  void   Put(int, double);

  /// Gets vector value
  double Get(int) const;

  //
  // Operators for vector access
  //

  /// Puts vector value
  double &operator()(int);

  /// Gets vector value
  double  operator()(int) const;

  //
  // Vector operators for doubles
  //

  /// Subtraction of a double
  irtkVector& operator-=(double);

  /// Addition of a double
  irtkVector& operator+=(double);

  /// Multiplication with a double
  irtkVector& operator*=(double);

  /// Division by a double
  irtkVector& operator/=(double);

  /// Return result of subtraction of a double
  irtkVector  operator- (double);

  /// Return result of addition of a double
  irtkVector  operator+ (double);

  /// Return result of multiplication with a double
  irtkVector  operator* (double);

  /// Return result of division by a double
  irtkVector  operator/ (double);

  //
  // Vector operators for vectors
  //

  /// Vector copy operator
  irtkVector& operator =(const irtkVector&);

  /// Vector subtraction operator
  irtkVector& operator-=(const irtkVector&);

  /// Vector addition operator
  irtkVector& operator+=(const irtkVector&);

  /// Vector componentwise multiplication operator (no scalar nor cross product)
  irtkVector& operator*=(const irtkVector&);

  /// Vector componentwise division operator
  irtkVector& operator/=(const irtkVector&);

  /// Return result for vector subtraction
  irtkVector  operator- (const irtkVector&);

  /// Return result for vector addition
  irtkVector  operator+ (const irtkVector&);

  /// Return result for componentwise vector multiplication (no scalar nor cross product)
  irtkVector  operator* (const irtkVector&);

  /// Return result for componentwise vector division
  irtkVector  operator/ (const irtkVector&);

  /// Comparison operator ==
  bool operator==(const irtkVector &);

#ifndef USE_STL
  /// Comparison operator != (if USE_STL is defined, negate == operator)
  bool operator!=(const irtkVector &);
#endif

  /// Comparison operator <
  bool operator<(const irtkVector &);

  /// Scalar/dot product
  double ScalarProduct(const irtkVector&);

  /// Vector/cross product
  irtkVector CrossProduct(const irtkVector&);

  /// Returns norm of a vector
  double Norm(void) const;

  /// Vector normalization
  void   Normalize(void);

  //
  // Vector in- and output
  //

  /// Interface to output stream
  friend ostream& operator<< (ostream&, const irtkVector&);

  /// Interface to input stream
  friend istream& operator>> (istream&, irtkVector&);

  /// Print vector
  void Print();

  /// Read vector from file
  void Read (char *);

  /// Write vector to file
  void Write(char *);


#ifdef USE_VXL
  /// Conversion to numerical recipes vector (memory must be allocated)
  template <class T> void Vector2Vnl(vnl_diag_matrix<T>*) const;

  /// Conversion from numerical recipes vector
  template <class T> void Vnl2Vector(vnl_diag_matrix<T>*);

#else
  /// Conversion to numerical recipes vector (memory must be allocated)
  void Vector2NR(float *) const;
  void Vector2NR(double *) const;

  /// Conversion from numerical recipes vector
  void NR2Vector(float *);
  void NR2Vector(double *);
#endif

};

//
// Access operators
//

inline int irtkVector::Rows() const
{
  return _rows;
}

inline void irtkVector::Put(int rows, double vector)
{
  _vector[rows] = vector;
}

inline double irtkVector::Get(int rows) const
{
  return _vector[rows];
}

inline double& irtkVector::operator()(int rows)
{
  return _vector[rows];
}

inline double irtkVector::operator()(int rows) const
{
  return _vector[rows];
}

//
// Vector operators for doubles
//

inline irtkVector& irtkVector::operator-=(double x)
{
  int i;

  for (i = 0; i < _rows; i++) {
    _vector[i] -= x;
  }

  return *this;
}

inline irtkVector& irtkVector::operator+=(double x)
{
  int i;

  for (i = 0; i < _rows; i++) {
    _vector[i] += x;
  }

  return *this;
}

inline irtkVector& irtkVector::operator*=(double x)
{
  int i;

  for (i = 0; i < _rows; i++) {
    _vector[i] *= x;
  }

  return *this;
}

inline irtkVector& irtkVector::operator/=(double x)
{
  int i;

  for (i = 0; i < _rows; i++) {
    _vector[i] /= x;
  }

  return *this;
}

inline irtkVector irtkVector::operator- (double x)
{
  int i;

  irtkVector m(_rows);
  for (i = 0; i < _rows; i++) {
    m._vector[i] = _vector[i] - x;
  }
  return m;
}

inline irtkVector irtkVector::operator+ (double x)
{
  int i;

  irtkVector m(_rows);
  for (i = 0; i < _rows; i++) {
    m._vector[i] = _vector[i] + x;
  }
  return m;
}

inline irtkVector irtkVector::operator* (double x)
{
  int i;

  irtkVector m(_rows);
  for (i = 0; i < _rows; i++) {
    m._vector[i] = _vector[i] * x;
  }
  return m;
}

inline irtkVector irtkVector::operator/ (double x)
{
  int i;

  irtkVector m(_rows);
  for (i = 0; i < _rows; i++) {
    m._vector[i] = _vector[i] / x;
  }
  return m;
}

//
// Vector operators for vectors
//

inline irtkVector& irtkVector::operator =(const irtkVector& v)
{
  int i;

  // Copy size
  this->Initialize(v._rows);

  // Copy vector
  for (i = 0; i < _rows; i++) {
    _vector[i] = v._vector[i];
  }
  return *this;
}

inline irtkVector& irtkVector::operator-=(const irtkVector& v)
{
  int i;

  if (_rows != v._rows) {
    cerr << "irtkVector::operator-=: Size mismatch" << endl;
    exit(1);
  }
  for (i = 0; i < _rows; i++) {
    _vector[i] -= v._vector[i];
  }

  return *this;
}

inline irtkVector& irtkVector::operator+=(const irtkVector& v)
{
  int i;

  if (_rows != v._rows) {
    cerr << "irtkVector::operator+=: Size mismatch" << endl;
    exit(1);
  }
  for (i = 0; i < _rows; i++) {
    _vector[i] += v._vector[i];
  }

  return *this;
}

inline irtkVector& irtkVector::operator*=(const irtkVector& v)
{
  int i;

  if (_rows != v._rows) {
    cerr << "irtkVector::operator*=: Size mismatch" << endl;
    exit(1);
  }
  for (i = 0; i < _rows; i++) {
    _vector[i] *= v._vector[i];
  }

  return *this;
}

inline irtkVector& irtkVector::operator/=(const irtkVector& v)
{
  int i;

  if (_rows != v._rows) {
    cerr << "irtkVector::operator/=: Size mismatch" << endl;
    exit(1);
  }
  for (i = 0; i < _rows; i++) {
    _vector[i] /= v._vector[i];
  }

  return *this;
}

inline irtkVector irtkVector::operator- (const irtkVector& v)
{
  int i;

  if (_rows != v._rows) {
    cerr << "irtkVector::operator-: Size mismatch" << endl;
    exit(1);
  }
  irtkVector m(_rows);
  for (i = 0; i < _rows; i++) {
    m._vector[i] = _vector[i] - v._vector[i];
  }
  return m;
}

inline irtkVector irtkVector::operator+ (const irtkVector& v)
{
  int i;

  if (_rows != v._rows) {
    cerr << "irtkVector::operator+: Size mismatch" << endl;
    exit(1);
  }
  irtkVector m(_rows);
  for (i = 0; i < _rows; i++) {
    m._vector[i] = _vector[i] + v._vector[i];
  }
  return m;
}

inline irtkVector irtkVector::operator* (const irtkVector& v)
{
  int i;

  if (_rows != v._rows) {
    cerr << "irtkVector::operator*: Size mismatch" << endl;
    exit(1);
  }
  irtkVector m(_rows);
  for (i = 0; i < _rows; i++) {
    m._vector[i] = _vector[i] * v._vector[i];
  }
  return m;
}

inline irtkVector irtkVector::operator/ (const irtkVector& v)
{
  int i;

  if (_rows != v._rows) {
    cerr << "irtkVector::operator/: Size mismatch" << endl;
    exit(1);
  }
  irtkVector m(_rows);
  for (i = 0; i < _rows; i++) {
    m._vector[i] = _vector[i] / v._vector[i];
  }
  return m;
}

//
// Comparison
//

inline bool  irtkVector::operator==(const irtkVector &v)
{
  if (_rows != v._rows) {
    return false;
  }
  for (int i = 0; i < _rows; i++) {
    if (_vector[i] != v._vector[i]) return false;
  }
  return true;

}

#ifndef USE_STL
inline bool  irtkVector::operator!=(const irtkVector &v)
{
  if (_rows != v._rows) {
    return true;
  }
  for (int i = 0; i < _rows; i++) {
    if (_vector[i] != v._vector[i]) return true;
  }
  return false;
}
#endif

inline bool  irtkVector::operator<(const irtkVector &v)
{
  if (_rows > v._rows) {
    return false;
  }
  for (int i = 0; i < _rows; i++) {
    if (_vector[i] >= v._vector[i]) return false;
  }
  return true;

}

//
// Vector products
//

inline double irtkVector::ScalarProduct(const irtkVector &v)
{
  int i;
  double scalar_product=0;

  if (_rows != v._rows) {
    cerr << "irtkVector::ScalarProduct: Size mismatch" << endl;
    exit(1);
  }

  for (i = 0; i < _rows; i++) {
    scalar_product += _vector[i]*v._vector[i];
  }
  return scalar_product;
}

inline irtkVector irtkVector::CrossProduct(const irtkVector &v)
{
  int i;

  if (_rows != v._rows) {
    cerr << "irtkVector::CrossProduct: Size mismatch" << endl;
    exit(1);
  }

  irtkVector m(_rows);
  for (i = 0; i < _rows; i++) {
    m._vector[i] = (_vector[(i+1)%_rows]*v._vector[(i+2)%_rows]-
                    _vector[(i+2)%_rows]*v._vector[(i+1)%_rows]);
  }
  return m;
}


//
// Functions
//

inline double irtkVector::Norm(void) const
{
  double norm = 0;

  for (int i = 0; i < _rows; i++) {
    norm += _vector[i]*_vector[i];
  }
  return sqrt(norm);
}

inline void irtkVector::Normalize(void)
{
  double norm = Norm();

  if (norm != 0) {
    *this /= norm;
  }
}

#ifdef USE_VXL

template <class T>
inline void irtkVector::Vector2Vnl(vnl_diag_matrix<T>* m) const
{
  unsigned i;

  for (i = 0; i < (unsigned) _rows; i++) {
    (*m)(i) = (T) _vector[i];
  }
}

template <class T>
inline void irtkVector::Vnl2Vector(vnl_diag_matrix<T>* m)
{
  unsigned i;

  for (i = 0; i < (unsigned) _rows; i++) {
    _vector[i] = (float) (*m)(i);
  }
}

#endif

#endif



