/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkTransformation.h>
#include <irtkRegistration.h>

#include <Fl_RView.h>

#include <Fl_RViewUI.h>

#include <Fl_RView_Histogram.h>

extern Fl_RViewUI *rviewUI;

Fl_HistogramWindow::Fl_HistogramWindow(int x, int y, int w, int h, const char *name, irtkRView  *viewer) : Fl_Window(x, y, w, h, name), _histogramWindow(viewer)
{
  _v = viewer;
}

Fl_HistogramWindow::~Fl_HistogramWindow()
{
}

void Fl_HistogramWindow::recalculate()
{
  _histogramWindow.CalculateHistograms();
}

void Fl_HistogramWindow::draw()
{
  int i, j;
  char buffer[256], buffer2[256];
  unsigned char r, g, b;
  double x, y, ox, oy;

  // Clear everything
  make_current();
  fl_draw_box(FL_FLAT_BOX, 0, 0, w(), h(), FL_WHITE);

  // Draw horizontal axis
  fl_color(225,225,225);
  fl_line_style(FL_DOT, 0);
  fl_line(10, round((h()-30)*0.00)+10, w()-10, round((h()-30)*0.00+10));
  fl_line(10, round((h()-30)*0.25)+10, w()-10, round((h()-30)*0.25+10));
  fl_line(10, round((h()-30)*0.50)+10, w()-10, round((h()-30)*0.50+10));
  fl_line(10, round((h()-30)*0.75)+10, w()-10, round((h()-30)*0.75+10));
  fl_line(10, round((h()-30)*1.00)+10, w()-10, round((h()-30)*1.00+10));

  // Draw vertical axis and labels
  for (i = 0; i <= HISTOGRAM_BINS; i += 32) {
    sprintf(buffer, "%d", i);
    fl_color(128,128,128);
    position(i, -1, x, y);
    fl_line(round(x), 10, round(x), h()-20);
    fl_draw(buffer, round(x)-5, h()-5);
    sprintf(buffer2, "%d", round(_histogramWindow._globalHistogram.BinToVal(i)));
    fl_draw(buffer2, round(x), h()-25);

  }

  // Compute maximum in histogram
  _maxHistogram = 0;
  for (int i = 0; i < HISTOGRAM_BINS; i++) {
    if (_histogramWindow._globalHistogram(i) > _maxHistogram) _maxHistogram = _histogramWindow._globalHistogram(i);
  }

  // Set up FL drawing colour and style
  fl_color(0,0,0);
  fl_line_style(FL_SOLID, 2);

  // Draw global histogram
  position(0,-1,x,y);
  for (int i = 1; i < HISTOGRAM_BINS; ++i) {
    ox = x;
    oy = y;
    position(i, -1, x, y);
    fl_line(round(ox), round(oy), round(x), round(y));
  }

  // Draw histogram for each structure
  for (j = 0; j < SHRT_MAX+1; j++) {

    // Check if structure is visible
    if ((_v->GetSegmentTable()->GetVisibility(j)) && (_v->GetSegmentTable()->IsValid(j))) {

      // Set up FL drawing colour style
      _v->GetSegmentTable()->GetColor(j,&r,&g,&b);
      fl_color(r,g,b);
      fl_line_style(FL_SOLID, 0);

      // Draw histogram
      position(0,j,x,y);
      for (int i = 1; i < HISTOGRAM_BINS; ++i) {
        ox = x;
        oy = y;
        position(i, j, x, y);
        fl_line(round(ox), round(oy), round(x), round(y));
      }
    }
  }
}

void Fl_HistogramWindow::position(int nbin, int nhistogram, double& x, double& y)
{
  x = nbin;
  if (nhistogram == -1) {
    y = _histogramWindow._globalHistogram(nbin);
  } else {
    y = _histogramWindow._localHistogram[nhistogram](nbin);
  }

  x = (x/HISTOGRAM_BINS)*(w() - 20) + 10;
  y = ((_maxHistogram - y) /_maxHistogram)*(h() - 30) + 10;
}

