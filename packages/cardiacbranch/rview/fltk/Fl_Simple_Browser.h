/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef FL_SIMPLE_BROWSER_H

#define FL_SIMPLE_BROWSER_H

#include <FL/Fl.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Hold_Browser.H>

/// This class provides a simple alternative to Fl_Browser
class Fl_Simple_Browser : public Fl_Hold_Browser
{

public:

  /// Normal FLTK constructor
  Fl_Simple_Browser(int x, int y, int w, int h, const char *label = 0);

  /// Override of Fl_Hold_Browser::handle()
  int handle(int);

};

inline Fl_Simple_Browser::Fl_Simple_Browser(int x, int y, int w, int h, const char *label) : Fl_Hold_Browser(x, y, w, h, label)
{
}

#endif
