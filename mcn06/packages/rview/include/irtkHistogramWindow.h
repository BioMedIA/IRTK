/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKHISTOGRAMWINDOW_H

#define _IRTKHISTOGRAMWINDOW_H

#include <irtkHistogram_1D.h>

#include <irtkRView.h>

#define HISTOGRAM_BINS 256

class irtkHistogramWindow
{

  friend class Fl_HistogramWindow;

protected:

  /// Pointer to registration viewer
  irtkRView *_v;

  /// Global histogram for entire image
  irtkHistogram_1D _globalHistogram;

  /// Global histogram for single segmentation
  irtkHistogram_1D _localHistogram[SHRT_MAX+1];

public:

  /// Constructor
  irtkHistogramWindow(irtkRView *);

  /// Destructor
  virtual ~irtkHistogramWindow() {};

  /// Compute histograms for everything
  void CalculateHistograms();

protected:

  /// Compute histogram for single segmentation
  void CalculateHistogram(int);
};

#endif
