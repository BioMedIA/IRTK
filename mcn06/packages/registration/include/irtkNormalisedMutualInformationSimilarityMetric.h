/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKNORMALISEDMUTUALINFORMATIONSIMILARITYMETRIC_H

#define _IRTKNORMALISEDMUTUALINFORMATIONSIMILARITYMETRIC_H

/**
 * Class for voxel similarity measure based on mutual information
 *
 */

class irtkNormalisedMutualInformationSimilarityMetric : public irtkHistogramSimilarityMetric
{

public:

  /// Constructor
  irtkNormalisedMutualInformationSimilarityMetric(int = 64, int = 64);

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkNormalisedMutualInformationSimilarityMetric::irtkNormalisedMutualInformationSimilarityMetric(int nbins_x, int nbins_y) : irtkHistogramSimilarityMetric (nbins_x, nbins_y)
{
}

inline double irtkNormalisedMutualInformationSimilarityMetric::Evaluate()
{
  return _histogram->NormalizedMutualInformation();
}

#endif
