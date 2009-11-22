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

ostream& operator<< (ostream& o, const irtkPoint &p)
{
  return o << p._x <<  " " << p._y << " " << p._z;
}

istream& operator>> (istream& i, irtkPoint &p)
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
  return i >> p._x >> p._y >> p._z;
}

