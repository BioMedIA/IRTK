/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKITERATIVERESAMPLING_H

#define _IRTKITERATIVERESAMPLING_H

#include <irtkResampling.h>

/**
 * Class for resampling of images (iterative)
 *
 * Iterative implementation of resampling for more accurate results
 *
 */

template <class VoxelType> class irtkIterativeResampling : public irtkResampling<VoxelType>
{

protected:

  irtkResampling<VoxelType> *_innerResample;

public:

  /// Constructor
  irtkIterativeResampling(double, double, double);

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Run the resampling filter
  virtual void   Run();

  virtual void   Run2D();

  virtual void   Run3D();
};

#endif
