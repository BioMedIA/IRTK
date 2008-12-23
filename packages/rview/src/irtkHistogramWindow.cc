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

#ifdef __APPLE__
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <irtkRView.h>

#ifdef HAS_VTK
#include <vtkPolyDataReader.h>
#include <vtkStructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#endif

#include <irtkHistogramWindow.h>

irtkHistogramWindow::irtkHistogramWindow(irtkRView  *viewer)
{
  _v = viewer;
}

void irtkHistogramWindow::CalculateHistograms()
{
  int i;

  if (_v->GetTarget()->IsEmpty()) {
    cerr << "No target image loaded." << endl;
    return;
  }

  // Calculate global histogram
  i = -1;
  CalculateHistogram(i);

  // Calculate histogram for each structure
  for (i = 0; i < _v->GetSegmentTable()->Size(); i++) {
    if ( _v->GetSegmentTable()->IsValid(i) == True) CalculateHistogram(i);
  }
}

void irtkHistogramWindow::CalculateHistogram(int label_id)
{
  int i;
  irtkGreyPixel min, max, *ptr2img, *ptr2seg;

  if (_v->GetTarget()->IsEmpty()) {
    cerr<< "No target image loaded." << endl;
    return;
  }

  _v->GetTarget()->GetMinMax(&min, &max);

  if (label_id < 0) {
    _globalHistogram.PutMin(min);
    _globalHistogram.PutMax(max);
    _globalHistogram.PutNumberOfBins(HISTOGRAM_BINS);

    ptr2img = _v->GetTarget()->GetPointerToVoxels();
    for (i = 0; i < _v->GetTarget()->GetNumberOfVoxels(); i++) {
      if (*ptr2img!=0) _globalHistogram.AddSample(*ptr2img);
      ptr2img++;
    }
  } else {
    if (_v->GetSegmentTable()->IsValid(label_id) == True) {
      _localHistogram[label_id].PutMin(min);
      _localHistogram[label_id].PutMax(max);
      _localHistogram[label_id].PutNumberOfBins(HISTOGRAM_BINS);
      ptr2img = _v->GetTarget()->GetPointerToVoxels();
      ptr2seg = _v->GetSegmentation()->GetPointerToVoxels();
      for (i = 0; i < _v->GetTarget()->GetNumberOfVoxels(); i++) {
        if ((*ptr2seg == label_id)&&(*ptr2img!=0)) _localHistogram[label_id].AddSample(*ptr2img);
        ptr2img++;
        ptr2seg++;
      }
    }
  }
}

