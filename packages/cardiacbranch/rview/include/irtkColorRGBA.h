/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCOLORRGBA_H

#define _IRTKCOLORRGBA_H

#define HAS_COLOR

#ifndef HAS_COLOR

typedef unsigned char irtkColorRGBA;

#else

class irtkColorRGBA
{

public:

  unsigned char r;
  unsigned char g;
  unsigned char b;
  float a;

  // Constructor (default)
  irtkColorRGBA();

  // Constructor (copy)
  irtkColorRGBA(const irtkColorRGBA &);

  // Copy operators
  irtkColorRGBA& operator= (const int &);

  // Copy operators
  irtkColorRGBA& operator= (const irtkColorRGBA &);

  // Convert HSV color to RGB color
  void HSVtoRGB(double, double, double);

  // Set RGB color
  void RGBtoHSV(double, double, double);

};

inline irtkColorRGBA::irtkColorRGBA()
{
  r = 0;
  g = 0;
  b = 0;
  a = 1;
}

inline irtkColorRGBA::irtkColorRGBA(const irtkColorRGBA &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
  a = c.a;
}

inline irtkColorRGBA& irtkColorRGBA::operator=(const irtkColorRGBA &c)
{
  r = c.r;
  g = c.g;
  b = c.b;
  a = c.a;
  return *this;
}

inline irtkColorRGBA& irtkColorRGBA::operator=(const int &c)
{
  r = c;
  g = c;
  b = c;
  return *this;
}

#endif

#endif
