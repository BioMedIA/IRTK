/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#ifndef _IRTKGRADIENTDESCENTOPTIMIZER_H

#define _IRTKGRADIENTDESCENTOPTIMIZER_H

/**
 * Generic class for gradient descent optimization of voxel-based
 * registration.
 */

class irtkGradientDescentOptimizer : public irtkOptimizer
{

public:

  /// Run the optimizer
  virtual double Run();

  /// Print name of the class
  virtual const char *NameOfClass();

};

inline const char *irtkGradientDescentOptimizer::NameOfClass()
{
  return "irtkGradientDescentOptimizer";
}

#endif
