/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSIMILARITYMETRIC_H

#define _IRTKSIMILARITYMETRIC_H

/**
 * Generic class for voxel similarity measures.
 *
 * This abstract class implements an interface to various voxel similarity
 * measures, including sums-of-squared differences, cross-correlation,
 * correlation ratio, mutual information and normalised mutual information.
 *
 */

class irtkSimilarityMetric : public irtkObject
{

public:

  /// Constructor
  irtkSimilarityMetric();

  /// Destructor
  virtual ~irtkSimilarityMetric();

  /// Copy metric
  static irtkSimilarityMetric *New(irtkSimilarityMetric *);

  /// Add sample
  virtual void Add(int, int, int =1) = 0;

  /// Remove sample
  virtual void Delete(int, int, int =1) = 0;

  /// Combine similarity metrics
  virtual void Combine(irtkSimilarityMetric *) = 0;

  /// Reset similarity metric
  virtual void Reset() = 0;

  /// Reset similarity metric
  virtual void Reset(irtkSimilarityMetric *) = 0;

  /// Evaluate similarity measure
  virtual double Evaluate() = 0;

};

inline irtkSimilarityMetric::irtkSimilarityMetric()
{
}

inline irtkSimilarityMetric::~irtkSimilarityMetric()
{
}

#include <irtkSSDSimilarityMetric.h>
#include <irtkCrossCorrelationSimilarityMetric.h>
#include <irtkMLSimilarityMetric.h>
#include <irtkHistogramSimilarityMetric.h>

inline irtkSimilarityMetric *irtkSimilarityMetric::New(irtkSimilarityMetric *metric)
{
  {
    irtkSSDSimilarityMetric *m = dynamic_cast<irtkSSDSimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkSSDSimilarityMetric;
    }
  }
  {
    irtkMutualInformationSimilarityMetric *m = dynamic_cast<irtkMutualInformationSimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkMutualInformationSimilarityMetric(m->NumberOfBinsX(), m->NumberOfBinsY());
    }
  }
  {
    irtkNormalisedMutualInformationSimilarityMetric *m = dynamic_cast<irtkNormalisedMutualInformationSimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkNormalisedMutualInformationSimilarityMetric(m->NumberOfBinsX(), m->NumberOfBinsY());
    }
  }
  {
    irtkJointEntropySimilarityMetric *m = dynamic_cast<irtkJointEntropySimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkJointEntropySimilarityMetric(m->NumberOfBinsX(), m->NumberOfBinsY());
    }
  }
  {
    irtkCrossCorrelationSimilarityMetric *m = dynamic_cast<irtkCrossCorrelationSimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkCrossCorrelationSimilarityMetric;
    }
  }
  {
    irtkKappaSimilarityMetric *m = dynamic_cast<irtkKappaSimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkKappaSimilarityMetric(m->NumberOfBinsX(), m->NumberOfBinsY());
    }
  }
  {
    irtkLabelConsistencySimilarityMetric *m = dynamic_cast<irtkLabelConsistencySimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkLabelConsistencySimilarityMetric;
    }
  }
  {
    irtkCorrelationRatioXYSimilarityMetric *m = dynamic_cast<irtkCorrelationRatioXYSimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkCorrelationRatioXYSimilarityMetric(m->NumberOfBinsX(), m->NumberOfBinsY());
    }
  }
  {
    irtkCorrelationRatioYXSimilarityMetric *m = dynamic_cast<irtkCorrelationRatioYXSimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkCorrelationRatioYXSimilarityMetric(m->NumberOfBinsX(), m->NumberOfBinsY());
    }
  }
  {
    irtkMLSimilarityMetric *m = dynamic_cast<irtkMLSimilarityMetric *>(metric);
    if (m != NULL) {
      return new irtkMLSimilarityMetric(m->_classification);
    }
  }

  cerr << "irtkSimilarityMetric::New: Failed" << endl;
  exit(1);
}


#endif
