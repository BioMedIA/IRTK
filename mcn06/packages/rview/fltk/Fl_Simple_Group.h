/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _FL_SIMPLE_GROUP_H

#define _FL_SIMPLE_GROUP_H

/* fltk includes */
#include <FL/Fl.H>
#include <FL/fl_draw.H>
#include <FL/Fl_Group.H>

//! This class provides a simple aesthetic alternative to Fl_Group
class Fl_Simple_Group : public Fl_Group
{

public:

  //! Normal FLTK constructor
  Fl_Simple_Group(int, int, int, int, const char *);

  //! Override of Fl_Group::draw()
  void draw();

};

#endif
