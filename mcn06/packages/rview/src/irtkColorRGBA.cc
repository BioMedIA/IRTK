/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkColorRGBA.h>

void irtkColorRGBA::RGBtoHSV(double cr, double cg, double cb)
{
  r = round(255*cr);
  g = round(255*cg);
  b = round(255*cb);
}

void irtkColorRGBA::HSVtoRGB(double h, double s, double v)
{
  int i;
  double f, p, q, t;

  if (s == 0) {
    if (h < 0) {
      r = round(255*v);
      g = round(255*v);
      b = round(255*v);
    } else {
      cerr << "irtkColorRGBA::HSV: Undefined HSV color" << endl;
      exit(1);
    }
  } else {
    if (h == 1) h = 0;
    h = h * 6;
    i = (int)floor(h);
    f = h - i;
    p = v * (1 - s);
    q = v * (1 - (s * f));
    t = v * (1 - (s * (1 - f)));
    switch (i) {
    case 0:
      r = round(255*v);
      g = round(255*t);
      b = round(255*p);
      break;
    case 1:
      r = round(255*q);
      g = round(255*v);
      b = round(255*p);
      break;
    case 2:
      r = round(255*p);
      g = round(255*v);
      b = round(255*t);
      break;
    case 3:
      r = round(255*p);
      g = round(255*q);
      b = round(255*v);
      break;
    case 4:
      r = round(255*t);
      g = round(255*p);
      b = round(255*v);
      break;
    case 5:
      r = round(255*v);
      g = round(255*p);
      b = round(255*q);
      break;
    default:
      break;
    }
  }
}

