/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <Fl_Simple_Group.h>

#include <string.h>
#include <stdio.h>

Fl_Simple_Group :: Fl_Simple_Group( int x, int y, int w, int h, const char *l )
    : Fl_Group( x, y, w, h, l )
{
  box( FL_EMBOSSED_FRAME );
  align( FL_ALIGN_LEFT | FL_ALIGN_INSIDE );
}

void Fl_Simple_Group :: draw()
{
  int lblW = 0, lblH, X;

  if ( label() == 0 )
    lblW = lblH = 0;
  else if ( strlen( label() ) == 0 )
    lblW = lblH = 0;
  else {
    measure_label( lblW, lblH );
    lblW += 4;
    lblH += 2;
  }

  // align the label
  if ( align() & FL_ALIGN_LEFT )
    X = 4;
  else if ( align() & FL_ALIGN_RIGHT )
    X = w() - lblW - 8;
  else
    X = w()/2 - lblW/2 - 2;

  // draw the main group box
  if ( damage() & ~FL_DAMAGE_CHILD )
    fl_draw_box( box(), x(), y()+lblH/2, w(), h()-lblH/2, color() );

  // clip and draw the children
  fl_clip( x()+2, y()+lblH+1, w()-4, h()-lblH-3 );
  draw_children();
  fl_pop_clip();

  // clear behind the label and draw it
  fl_color( color() );
  fl_rectf( x()+X, y(), lblW+4, lblH );
  fl_color( labelcolor() );
  draw_label( x()+X+2, y(), lblW, lblH, FL_ALIGN_CENTER );
}
