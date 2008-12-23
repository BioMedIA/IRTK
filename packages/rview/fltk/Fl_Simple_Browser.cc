/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <Fl_Simple_Browser.h>

#include <string.h>
#include <stdio.h>

int Fl_Simple_Browser::handle(int event)
{
  switch (event) {
  case FL_KEYBOARD:
  case FL_SHORTCUT:
    if (Fl::event_key() == FL_BackSpace) {
      do_callback();
      return 1;
    }
    return this->Fl_Browser_::handle(event);
  default:
    return this->Fl_Browser_::handle(event);
  }
  return 1;
}

