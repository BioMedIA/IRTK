/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _FL_RVIEW_H

#define _FL_RVIEW_H

#ifdef __APPLE__
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>

#include <irtkRView.h>


class Fl_RView : public Fl_Gl_Window
{

public:

  /// Pointer to the registration viewer
  irtkRView *v;

  /// Constructor
  Fl_RView(int, int, int, int, const char *);

  /// Default draw function
  void draw();

  /// Default function to handle events
  int  handle(int);

};

#endif
