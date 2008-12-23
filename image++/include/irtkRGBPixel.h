/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKRGBPIXEL_H

#define _IRTKRGBPIXEL_H

class irtkRGBPixel
{

public:

  /// Red
  unsigned char _r;

  /// Green
  unsigned char _g;

  /// Blue
  unsigned char _b;

  //
  // Constructors and destructor
  //

  /// Constructor
  irtkRGBPixel();

  /// Constructor with three coordinates
  irtkRGBPixel(unsigned char, unsigned char, unsigned char);

  /// Constructor with Point
  irtkRGBPixel(const irtkRGBPixel &);

  /// Default destructor
  virtual ~irtkRGBPixel(void);

  //
  // Operators for Point
  //

  /// Copy operator for point
  irtkRGBPixel& operator =(const irtkRGBPixel&);

  /// Return result of point substraction
  irtkRGBPixel  operator- (const irtkRGBPixel&);

  /// Return result of point addition
  irtkRGBPixel  operator+ (const irtkRGBPixel&);

  /// Return result of point multiplication
  irtkRGBPixel  operator* (const irtkRGBPixel&);

  /// Return result of point division
  irtkRGBPixel  operator/ (const irtkRGBPixel&);

  /// Substraction operator for point
  irtkRGBPixel& operator-=(const irtkRGBPixel&);

  /// Addition operator for point
  irtkRGBPixel& operator+=(const irtkRGBPixel&);

  /// Multiplication operator for point
  irtkRGBPixel& operator*=(const irtkRGBPixel&);

  /// Division operator for point
  irtkRGBPixel& operator/=(const irtkRGBPixel&);

  /// Comparison operator ==
  int    operator==(const irtkRGBPixel&);

  /// Comparison operator != (if USE_STL is defined, negate == operator)
  int    operator!=(const irtkRGBPixel&);

  /// Comparison operator <
  int    operator<(const irtkRGBPixel&);

  /// Comparison operator >
  int    operator>(const irtkRGBPixel&);

  /// Comparison operator >=
  int    operator>=(const irtkRGBPixel&);

};

inline irtkRGBPixel::irtkRGBPixel()
{
  _r = 0;
  _g = 0;
  _b = 0;
}

inline irtkRGBPixel::irtkRGBPixel(unsigned char x, unsigned char y, unsigned char z)
{
  _r = x;
  _g = y;
  _b = z;
}

inline irtkRGBPixel::irtkRGBPixel(const irtkRGBPixel& p)
{
  _r = p._r;
  _g = p._g;
  _b = p._b;
}

inline irtkRGBPixel::~irtkRGBPixel()
{
}

inline irtkRGBPixel& irtkRGBPixel::operator=(const irtkRGBPixel& p)
{
  _r = p._r;
  _g = p._g;
  _b = p._b;
  return *this;
}

inline irtkRGBPixel& irtkRGBPixel::operator+=(const irtkRGBPixel& p)
{
  _r += p._r;
  _g += p._g;
  _b += p._b;
  return *this;
}

inline irtkRGBPixel& irtkRGBPixel::operator-=(const irtkRGBPixel& p)
{
  _r -= p._r;
  _g -= p._g;
  _b -= p._b;
  return *this;
}

inline irtkRGBPixel& irtkRGBPixel::operator*=(const irtkRGBPixel& p)
{
  _r *= p._r;
  _g *= p._g;
  _b *= p._b;
  return *this;
}

inline irtkRGBPixel& irtkRGBPixel::operator/=(const irtkRGBPixel& p)
{
  _r /= p._r;
  _g /= p._g;
  _b /= p._b;
  return *this;
}

inline irtkRGBPixel irtkRGBPixel::operator+ (const irtkRGBPixel& p)
{
  irtkRGBPixel tmp;

  tmp._r = _r + p._r;
  tmp._g = _g + p._g;
  tmp._b = _b + p._b;
  return tmp;
}

inline irtkRGBPixel irtkRGBPixel::operator- (const irtkRGBPixel& p)
{
  irtkRGBPixel tmp;

  tmp._r = _r - p._r;
  tmp._g = _g - p._g;
  tmp._b = _b - p._b;
  return tmp;
}

inline irtkRGBPixel irtkRGBPixel::operator* (const irtkRGBPixel& p)
{
  irtkRGBPixel tmp;

  tmp._r = _r * p._r;
  tmp._g = _g * p._g;
  tmp._b = _b * p._b;
  return tmp;
}

inline irtkRGBPixel irtkRGBPixel::operator/ (const irtkRGBPixel& p)
{
  irtkRGBPixel tmp;

  tmp._r = _r / p._r;
  tmp._g = _g / p._g;
  tmp._b = _b / p._b;
  return tmp;
}

inline int irtkRGBPixel::operator==(irtkRGBPixel const &p)
{
  return ((_r == p._r) && (_g == p._g) && (_b == p._b));
}

inline int irtkRGBPixel::operator!=(irtkRGBPixel const &p)
{
  return ((_r != p._r) || (_g != p._g) || (_b != p._b));
}

inline int irtkRGBPixel::operator>(irtkRGBPixel const &p)
{
  return ((_r > p._r) && (_g > p._g) && (_b > p._b));
}

inline int irtkRGBPixel::operator<(irtkRGBPixel const &p)
{
  return ((_r < p._r) && (_g < p._g) && (_b < p._b));
}

inline int irtkRGBPixel::operator>=(irtkRGBPixel const &p)
{
  return ((_r >= p._r) && (_g >= p._g) && (_b >= p._b));
}
#endif
