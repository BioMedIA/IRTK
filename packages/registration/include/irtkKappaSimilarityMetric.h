/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKKAPPASIMILARITYMETRIC_H

#define _IRTKKAPPASIMILARITYMETRIC_H

/**
 * Class for voxel similarity measure based on the kappa measure
 *
 */

class irtkKappaSimilarityMetric : public irtkHistogramSimilarityMetric
{

public:

  /// Constructor
  irtkKappaSimilarityMetric(int = 64, int = 64);

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkKappaSimilarityMetric::irtkKappaSimilarityMetric(int nbins_x, int nbins_y) : irtkHistogramSimilarityMetric (nbins_x, nbins_y)
{
}

inline double irtkKappaSimilarityMetric::Evaluate()
{
  return _histogram->Kappa();
}

#endif
