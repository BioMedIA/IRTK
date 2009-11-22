/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCORRELATIONRATIOYXSIMILARITYMETRIC_H

#define _IRTKCORRELATIONRATIOYXSIMILARITYMETRIC_H

/**
 * Class for voxel similarity measure based on correlation ratio
 *
 */

class irtkCorrelationRatioYXSimilarityMetric : public irtkHistogramSimilarityMetric
{

public:

  /// Constructor
  irtkCorrelationRatioYXSimilarityMetric(int = 64, int = 64);

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkCorrelationRatioYXSimilarityMetric::irtkCorrelationRatioYXSimilarityMetric(int nbins_x, int nbins_y) : irtkHistogramSimilarityMetric (nbins_x, nbins_y)
{
}

inline double irtkCorrelationRatioYXSimilarityMetric::Evaluate()
{
  return _histogram->CorrelationRatioYX();
}

#endif
