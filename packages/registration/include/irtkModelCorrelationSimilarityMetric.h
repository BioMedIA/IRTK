/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMODELCORRELATIONSIMILARITYMETRIC_H

#define _IRTKMODELCORRELATIONSIMILARITYMETRIC_H

#ifdef HAS_VTK

/**
 * Generic class for model similarity measures.
 *
 */

class irtkModelCorrelationSimilarityMetric : public irtkModelSimilarityMetric
{

  /// Interpolator
  irtkInterpolateImageFunction *_interpolator;

  /// Domain which can be interpolated fast
  double _x1, _y1, _z1, _x2, _y2, _z2;

  /// Mean
  double _mean;

  /// Variance
  double _var;

  /// Number of samples in profile
  int _n;

  /// Step size for profile
  double _ds;

  /// Profile
  double *_profile;

public:

  /// Constructor
  irtkModelCorrelationSimilarityMetric(irtkGreyImage *, int n = 10, double ds = 1);

  /// Destructor
  virtual ~irtkModelCorrelationSimilarityMetric();

  /// Evaluate metric
  virtual double Evaluate(double *, double * = NULL, double * = NULL);

};

#define MAX_PROFILE 100

#endif

#endif
