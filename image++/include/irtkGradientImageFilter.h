/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGRADIENTIMAGEFILTER_H

#define _IRTKGRADIENTIMAGEFILTER_H

/**
 * Class for calculating the gradient of an image.
 *
 * The class provides an interface to calculating the gradient in the
 * x- , y- and z- directions.
 */

#include <irtkImageToImage.h>

template <class VoxelType> class irtkGradientImageFilter : public irtkImageToImage<VoxelType>
{

protected:

	/// Type of gradient
	int _type;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /** Returns whether the filter requires buffering. This filter requires
   *  buffering and returns 0.
   */
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  // Type of gradient vector to compute
  const static int GRADIENT_X          = 0;
  const static int GRADIENT_Y          = 1;
  const static int GRADIENT_Z          = 2;
  const static int GRADIENT_MAGNITUDE  = 3;
  const static int GRADIENT_VECTOR     = 4;

  /// Constructor
  irtkGradientImageFilter(int type = irtkGradientImageFilter::GRADIENT_MAGNITUDE);

  /// Run the convolution filter
  virtual void Run();
};

#endif
