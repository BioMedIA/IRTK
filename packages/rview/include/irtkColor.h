/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCOLOR_H

#define _IRTKCOLOR_H

#define HAS_COLOR

#ifndef HAS_COLOR

typedef unsigned char irtkColor;
typedef unsigned char irtkColorRGBA;

#else

class irtkColor
{

public:

  unsigned char r;
  unsigned char g;
  unsigned char b;

  // Constructor (default)
  irtkColor();

  // Constructor (copy)
  irtkColor(const irtkColor &);

  // Constructor (copy)
  irtkColor(const irtkColorRGBA &);

  // Copy operators
  irtkColor& operator= (const int &);

  // Copy operators
  irtkColor& operator= (const irtkColor &);

  // Copy operators
  irtkColor& operator= (const irtkColorRGBA &);

  // Convert HSV color to RGB color
  void HSVtoRGB(double, double, double);

  // Set RGB color
  void RGBtoHSV(double, double, double);

};

inline irtkColor::irtkColor()
{
  r = 0;
  g = 0;
  b = 0;
}

inline irtkColor::irtkColor(const irtkColor &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
}

inline irtkColor::irtkColor(const irtkColorRGBA &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
}

inline irtkColor& irtkColor::operator=(const irtkColor &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
  return *this;
}

inline irtkColor& irtkColor::operator=(const irtkColorRGBA &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
  return *this;
}

inline irtkColor& irtkColor::operator=(const int &c)
{
  r = c;
  g = c;
  b = c;
  return *this;
}

#endif

#endif
