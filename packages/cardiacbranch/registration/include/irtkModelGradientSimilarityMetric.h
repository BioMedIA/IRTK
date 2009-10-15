/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMODELGRADIENTSIMILARITYMETRIC_H

#define _IRTKMODELGRADIENTSIMILARITYMETRIC_H

#ifdef HAS_VTK

/**
 * Generic class for model similarity measures.
 *
 */

class irtkModelGradientSimilarityMetric : public irtkModelSimilarityMetric
{

	/// Gradient image (t = 0 contains the gradient magnitude, t = 1, ... 3 contains the gradient vector itself
	irtkRealImage *_gradient;
	
  /// Interpolator
  irtkInterpolateImageFunction *_interpolator;

  /// Domain which can be interpolated fast
  double _x1, _y1, _z1, _x2, _y2, _z2;

public:

  /// Constructor
  irtkModelGradientSimilarityMetric(irtkGreyImage *);

  /// Destructor
  virtual ~irtkModelGradientSimilarityMetric();

  /// Evaluate metric
  virtual double Evaluate(double *, double * = NULL, double * = NULL);

};

#endif

#endif
