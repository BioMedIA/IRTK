/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKHISTOGRAMSIMILARITYMETRIC_H

#define _IRTKHISTOGRAMSIMILARITYMETRIC_H

/**
 * Generic class for histogram-based voxel similarity measures.
 *
 * This class implements various voxel similarity measures, including
 * correlation ratio, mutual information and normalised mutual information.
 *
 */

class irtkHistogramSimilarityMetric : public irtkSimilarityMetric
{

protected:

  /// Histogram
  irtkHistogram_2D *_histogram;

public:

  /// Constructor
  irtkHistogramSimilarityMetric(int = 64, int = 64);

  /// Destructor
  ~irtkHistogramSimilarityMetric();

  /// Add sample
  virtual void Add(int, int);

  /// Remove sample
  virtual void Delete(int, int);

  /// Combine similarity metrics
  virtual void Combine(irtkSimilarityMetric *);

  /// Reset similarity metric
  virtual void Reset();

  /// Reset similarity metric
  virtual void Reset(irtkSimilarityMetric *);

  /// Return number of bins in X
  int NumberOfBinsX();

  /// Return number of bins in Y
  int NumberOfBinsY();

};

inline irtkHistogramSimilarityMetric::irtkHistogramSimilarityMetric(int nbins_x, int nbins_y)
{
  _histogram = new irtkHistogram_2D(nbins_x, nbins_y);
}

inline irtkHistogramSimilarityMetric::~irtkHistogramSimilarityMetric()
{
  delete _histogram;
}

inline void irtkHistogramSimilarityMetric::Add(int x, int y)
{
  _histogram->Add(x, y);
}

inline void irtkHistogramSimilarityMetric::Delete(int x, int y)
{
  _histogram->Delete(x, y);
}

inline void irtkHistogramSimilarityMetric::Combine(irtkSimilarityMetric *metric)
{
  int i, j;
  irtkHistogramSimilarityMetric *m = dynamic_cast<irtkHistogramSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkHistogramSimilarityMetric::Combine: Dynamic cast failed" << endl;
    exit(1);
  }

  if ((_histogram->GetNumberOfBinsX() != m->_histogram->GetNumberOfBinsX()) ||
      (_histogram->GetNumberOfBinsY() != m->_histogram->GetNumberOfBinsY())) {
    cerr << "irtkHistogramSimilarityMetric::Combine: Number of bins differs" << endl;
    exit(1);
  }
  for (j = 0; j < _histogram->GetNumberOfBinsY(); j++) {
    for (i = 0; i < _histogram->GetNumberOfBinsX(); i++) {
      _histogram->Add(i, j, m->_histogram->irtkHistogram_2D::operator()(i, j));
    }
  }
}

inline void irtkHistogramSimilarityMetric::Reset()
{
  _histogram->Reset();
}

inline void irtkHistogramSimilarityMetric::Reset(irtkSimilarityMetric *metric)
{
  irtkHistogramSimilarityMetric *m = dynamic_cast<irtkHistogramSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkHistogramSimilarityMetric::Reset: Dynamic cast failed" << endl;
    exit(1);
  }

  _histogram->Reset(*m->_histogram);
}

inline int irtkHistogramSimilarityMetric::NumberOfBinsX()
{
  return _histogram->GetNumberOfBinsX();
}

inline int irtkHistogramSimilarityMetric::NumberOfBinsY()
{
  return _histogram->GetNumberOfBinsY();
}

#include <irtkGenericHistogramSimilarityMetric.h>
#include <irtkMutualInformationSimilarityMetric.h>
#include <irtkNormalisedMutualInformationSimilarityMetric.h>
#include <irtkJointEntropySimilarityMetric.h>
#include <irtkCorrelationRatioXYSimilarityMetric.h>
#include <irtkCorrelationRatioYXSimilarityMetric.h>
#include <irtkKappaSimilarityMetric.h>
#include <irtkLabelConsistencySimilarityMetric.h>

#endif
