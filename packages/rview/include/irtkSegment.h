/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSEGMENT_H

#define _IRTKSEGMENT_H

#define HEX_LENGTH 7

#include <irtkImage.h>

#include <irtkColorRGBA.h>
#include <irtkColor.h>

class irtkSegment
{

  friend class irtkRView;
  friend class irtkViewer;

protected:

  /// Color
  irtkColor _color;

  /// Color in hex format
  char _hexColor[HEX_LENGTH];

  /// Name of structure
  char *_label;

  /// Transperancy
  double _trans;

  /// Visibility flag
  int _visible;

  void setHexColor(void);

public:

  // Constructor (default)
  irtkSegment();

  // Constructor (existing)
  irtkSegment(char*, unsigned char, unsigned char, unsigned char, double, int = true);

  // Destructor
  virtual ~irtkSegment(void);

  // Copy operator
  irtkSegment& operator = (const irtkSegment&);

  /// Set color
  void setColor(unsigned char, unsigned char, unsigned char);

  /// Set label
  void setLabel(char*);

  /// Set transperancy
  void setTrans(double);

  /// Set to visible
  void setVisibility(int);

  /// Return color
  void getColor(unsigned char*, unsigned char*, unsigned char*) const;

  /// Return color in hex
  void getHex(char *) const;

  /// Return transperancy
  double getTrans() const;

  /// Return label
  char *getLabel() const;

  /// Return if visible
  int   getVisibility() const;

  // General Methods
  void  rgb2Hex(int, int, int, char*);
  char* int2Hex(int, unsigned char);

};

#endif
