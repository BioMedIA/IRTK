/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkSegment.h>

irtkSegment::irtkSegment()
{
  // Visibility
  _visible = False;

  // Label name
  _label = NULL;

  // Colour
  setColor(0, 0, 0);

  // Transperancy
  _trans = 0.0;
}

irtkSegment::irtkSegment(char* label, unsigned char cr, unsigned char cg, unsigned char cb, double t, int v)
{
  // Visibility
  _visible = v;

  // Label name
  _label = strdup(label);

  // Colour
  setColor(cr, cg, cb);

  // Transperancy
  _trans = t;
}

irtkSegment::~irtkSegment()
{
  if (_label != NULL) free(_label);
}

irtkSegment& irtkSegment::operator =(const irtkSegment& s)
{
  strncpy(_hexColor, s._hexColor, HEX_LENGTH);
  _color = s._color;
  _label = strdup(s._label);
  _trans = s._trans;
  return *this;
}

void irtkSegment::setLabel(char* label)
{
  if (_label != NULL) free(_label);
  if (label != NULL) {
    _label = strdup(label);
  } else {
    _label = NULL;
  }
}

void irtkSegment::setColor(unsigned char r, unsigned char g, unsigned char b)
{
  _color.r = r;
  _color.g = g;
  _color.b = b;
  this->setHexColor();
}

void irtkSegment::setHexColor(void)
{
  this->rgb2Hex(_color.r, _color.g, _color.b, _hexColor);
}

void irtkSegment::setTrans(double t)
{
  _trans = t;
}

void irtkSegment::setVisibility(int vis)
{
  _visible = vis;
}

void irtkSegment::getColor(unsigned char *r, unsigned char *g, unsigned char *b) const
{
  *r = _color.r;
  *g = _color.g;
  *b = _color.b;
}

void irtkSegment::getHex(char *hexColor) const
{
  strncpy(hexColor, _hexColor, HEX_LENGTH);
}

double irtkSegment::getTrans() const
{
  return _trans;
}

char *irtkSegment::getLabel() const
{
  return _label;
}

int irtkSegment::getVisibility() const
{
  return _visible;
}


void irtkSegment::rgb2Hex(int r, int g, int b, char* h)
{
  char* hr = int2Hex(r, 2);
  char* hg = int2Hex(g, 2);
  char* hb = int2Hex(b, 2);

  h[0] = hr[1];
  h[1] = hr[2];
  h[2] = hg[1];
  h[3] = hg[2];
  h[4] = hb[1];
  h[5] = hb[2];
  h[6] = '\0';

  delete [] hr;
  delete [] hg;
  delete [] hb;
}

char* irtkSegment::int2Hex(int n, unsigned char round)
{
  unsigned char size = round;
  char i, *hex, *hexLookup = "0123456789ABCDEF";
  int temp = n;

  hex = new char[size+2];

  if (n<0) {
    hex[0]='-';
    n *= -1;
  } else {
    hex[0]=' ';
  }

  char mask = 0x000f;

  for (i=0; i<size; i++) {
    temp = n;
    temp >>=(4*i);
    temp &= mask;
    hex[size-i]=hexLookup[temp];
  }

  hex[size+1]= '\0';

  return hex;
}
