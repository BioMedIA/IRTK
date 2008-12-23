/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSTEEPESTGRADIENTDESCENTOPTIMIZER_H

#define _IRTKSTEEPESTGRADIENTDESCENTOPTIMIZER_H

/**
 * Generic class for steepest gradient descent optimization of voxel-based
 * registration.
 */

class irtkSteepestGradientDescentOptimizer : public irtkOptimizer
{

public:

  /// Run the optimizer
  virtual double Run();

  /// Print name of the class
  virtual const char *NameOfClass();

};

inline const char *irtkSteepestGradientDescentOptimizer::NameOfClass()
{
  return "irtkSteepestGradientDescentOptimizer";
}

#endif
