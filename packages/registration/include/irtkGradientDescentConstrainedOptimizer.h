/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGRADIENTDESCENTCONSTRAINEDOPTIMIZER_H

#define _IRTKGRADIENTDESCENTCONSTRAINEDOPTIMIZER_H

/**
 * Class for gradient descent optimization of voxel-based
 * registration which enforces hard limits on the registration
 * parameters
 */

class irtkGradientDescentConstrainedOptimizer : public irtkOptimizer
{

  /// Hard limits
  double _limits;

public:

  /// Constructor
  irtkGradientDescentConstrainedOptimizer();

  /// Run the optimizer
  virtual double Run();

  /// Print name of the class
  virtual const char *NameOfClass();

  /// Set hard limits for parameters
  virtual void SetLimits(double);

};

inline const char *irtkGradientDescentConstrainedOptimizer::NameOfClass()
{
  return "irtkGradientDescentOptimizer";
}

inline void irtkGradientDescentConstrainedOptimizer::SetLimits(double limits)
{
  _limits = limits;
}

#endif
