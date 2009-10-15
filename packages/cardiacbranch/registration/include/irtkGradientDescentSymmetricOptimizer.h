/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGRADIENTDESCENTSYMMETRICOPTIMIZER_H

#define _IRTKGRADIENTDESCENTSYMMETRICOPTIMIZER_H

/**
 * Generic class for gradient descent optimization of voxel-based
 * symmetric registration.
 */

class irtkGradientDescentSymmetricOptimizer : public irtkSymmetricOptimizer
{

public:

  /// Run the optimizer
  virtual double Run();

  /// Print name of the class
  virtual const char *NameOfClass();

};

inline const char *irtkGradientDescentSymmetricOptimizer::NameOfClass()
{
  return "irtkGradientDescentSymmetricOptimizer";
}

#endif
