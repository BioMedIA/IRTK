/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKJOINTENTROPYSIMILARITYMETRIC_H

#define _IRTKJOINTENTROPYSIMILARITYMETRIC_H

/**
 * Class for voxel similarity measure based on joint entropy
 *
 */

class irtkJointEntropySimilarityMetric : public irtkHistogramSimilarityMetric
{

public:

  /// Constructor
  irtkJointEntropySimilarityMetric(int = 64, int = 64);

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkJointEntropySimilarityMetric::irtkJointEntropySimilarityMetric(int nbins_x, int nbins_y) : irtkHistogramSimilarityMetric (nbins_x, nbins_y)
{
}

inline double irtkJointEntropySimilarityMetric::Evaluate()
{
  return -_histogram->JointEntropy();
}

#endif
