/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGRADIENTIMAGE_H

#define _IRTKGRADIENTIMAGE_H

/**
 * Class for caluclating the gradient of an image.
 * The class provides an iterface to subclasses which calculate the gradient in
 * x- , y- and z- directions.
 */

#include <irtkImageToImage.h>

template <class VoxelType> class irtkGradientImage : public irtkImageToImage<VoxelType>
{

protected:

  /** Returns whether the filter requires buffering. This filter requires
   *  buffering and returns 0.
   */
  virtual Bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  // Runs the filter on a single voxel.
  virtual double Run(int, int, int, int);

public:

  /// Run the convolution filter
  virtual void Run();
};

#include <irtkGradientImageX.h>
#include <irtkGradientImageY.h>
#include <irtkGradientImageZ.h>

#endif
