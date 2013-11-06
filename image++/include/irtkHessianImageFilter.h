/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKHESSIANIMAGEFILTER_H

#define _IRTKHESSIANIMAGEFILTER_H

/**
 * Class for calculating the second order gradient of an image.
 *
 * The class provides an interface to calculating the second order gradient in the
 * x-x-, x-y-, x-z-, y-y-, y-z- and z-z- directions.
 */

#include <irtkImageToImage.h>

template <class VoxelType> class irtkHessianImageFilter : public irtkImageToImage<VoxelType>
{

protected:

  /// Type of gradient
  int _type;

  /// Padding
  int _Padding;

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
  const static int HESSIAN_XX = 0;
  const static int HESSIAN_XY = 1;
  const static int HESSIAN_XZ = 2;
  const static int HESSIAN_YY = 3;
  const static int HESSIAN_YZ = 4;
  const static int HESSIAN_ZZ = 5;
  const static int HESSIAN_VECTOR = 6;

  /// Constructor
  irtkHessianImageFilter(int type = irtkHessianImageFilter::HESSIAN_VECTOR);

  /// Run the convolution filter
  virtual void Run();

  /// Set Padding
  virtual SetMacro(Padding,VoxelType);
};

#endif
