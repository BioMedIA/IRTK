/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCORRELATIONRATIOXYSIMILARITYMETRIC_H

#define _IRTKCORRELATIONRATIOXYSIMILARITYMETRIC_H

/**
 * Class for voxel similarity measure based on correlation ratio
 *
 */

class irtkCorrelationRatioXYSimilarityMetric : public irtkHistogramSimilarityMetric
{

public:

  /// Constructor
  irtkCorrelationRatioXYSimilarityMetric(int = 64, int = 64);

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkCorrelationRatioXYSimilarityMetric::irtkCorrelationRatioXYSimilarityMetric(int nbins_x, int nbins_y) : irtkHistogramSimilarityMetric (nbins_x, nbins_y)
{
}

inline double irtkCorrelationRatioXYSimilarityMetric::Evaluate()
{
  return _histogram->CorrelationRatioXY();
}

#endif
