/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMODELSIMILARITYMETRIC_H

#define _IRTKMODELSIMILARITYMETRIC_H

#ifdef HAS_VTK

/**
 * Generic class for model similarity measures.
 *
 */

class irtkModelSimilarityMetric : public irtkObject
{

protected:
	
	/// Target image
	irtkGreyImage *_image;
		
public:

  /// Constructor
  irtkModelSimilarityMetric(irtkGreyImage *);

  /// Destructor
  virtual ~irtkModelSimilarityMetric();

  /// Evaluate metric
  virtual double Evaluate(double *, double * = NULL, double * = NULL) = 0;

};

inline irtkModelSimilarityMetric::irtkModelSimilarityMetric(irtkGreyImage *image)
{
	_image  = image;
}

inline irtkModelSimilarityMetric::~irtkModelSimilarityMetric()
{
}

#include <irtkModelGradientSimilarityMetric.h>
#include <irtkModelCorrelationSimilarityMetric.h>

#endif

#endif
