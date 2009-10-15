/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKDOWNHILLDESCENTOPTIMIZER_H

#define _IRTKDOWNHILLDESCENTOPTIMIZER_H

/**
 * Generic class for downhill descent optimization of voxel-based
 * registration.
 */

class irtkDownhillDescentOptimizer : public irtkOptimizer
{

public:

  /// Optimization method
  virtual double Run();

  /// Print name of the class
  virtual const char *NameOfClass();
};

inline const char *irtkDownhillDescentOptimizer::NameOfClass()
{
  return "irtkDownhillGradientDescentOptimizer";
}

#endif
