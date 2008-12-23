/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _FL_HISTOGRAMWINDOW_H

#define _FL_HISTOGRAMWINDOW_H

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/fl_draw.H>

#include <irtkRView.h>

class Fl_HistogramWindow : public Fl_Window
{

protected:

  /// Pointer to histogram window
  irtkHistogramWindow _histogramWindow;

  /// Maximum in histogram
  int _maxHistogram;

public:

  /// Pointer to the registration viewer
  irtkRView *_v;

  /// Constructor
  Fl_HistogramWindow(int, int, int, int, const char *, irtkRView *);

  /// Destructor
  ~Fl_HistogramWindow();

  /// Default draw function
  void draw();

  /// Compute position
  void position(int, int, double& , double&);

  /// Recalculate histogram
  void recalculate();

};

#endif
