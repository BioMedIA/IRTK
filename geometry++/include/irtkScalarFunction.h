/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSCALARFUNCTION_H

#define _IRTKSCALARFUNCTION_H

/**

  Scalar function class.

*/

class irtkScalarFunction : public irtkObject
{

public:

  /// Virtual destructor
  virtual ~irtkScalarFunction();

  /// Evaluation function (pure virtual)
  virtual double Evaluate(double x, double y, double z) = 0;
};

#include <irtkScalarGaussian.h>

#endif
