/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkGeometry.h>

#include <irtkPoint4D.h>

ostream& operator<< (ostream& o, const irtkPoint4D &p)
{
  return o << p._x <<  " " << p._y << " " << p._z<<p._t;
}

istream& operator>> (istream& i, irtkPoint4D &p)
{
  char c;

  c = i.peek();
  while ((c == ',') || (c == ' ') || (c == '\t') || (c == '\n')) {
    c = i.get();
    if (!i.eof()) {
      c = i.peek();
    } else {
      break;
    }
  }
  return i >> p._x >> p._y >> p._z >> p._t;
}

