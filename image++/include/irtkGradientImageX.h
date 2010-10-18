/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGRADIENTIMAGE_X_H

#define _IRTKGRADIENTIMAGE_X_H

#include <irtkGradientImage.h>

/**
 * Class for caluclating the gradient of an image in the x-direction.
 */

template <class VoxelType> class irtkGradientImageX : public irtkGradientImage<VoxelType>
{

protected:


  /** Returns whether the filter requires buffering. This filter requires
   *  buffering and returns 0.
   */
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  // Calculate the gradient on a single voxel.
  virtual double Run(int, int, int, int);

public:

  /// Run the convolution filter
  virtual void Run();
};


#endif
