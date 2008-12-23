/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGAUSSIANBLURRING_H

#define _IRTKGAUSSIANBLURRING_H

#include <irtkImageToImage.h>

/**
 * Class for Gaussian blurring of images
 *
 * This class defines and implements the Gaussian blurring of images. The
 * blurring is implemented by three successive 1D convolutions with a 1D
 * Gaussian kernel.
 */

template <class VoxelType> class irtkGaussianBlurring : public irtkImageToImage<VoxelType>
{

protected:

  /// Sigma (standard deviation of Gaussian kernel)
  double _Sigma;

  /// Returns whether the filter requires buffering
  virtual Bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor
  irtkGaussianBlurring(double);

  /// Destructor
  ~irtkGaussianBlurring();

  /// Run Gaussian blurring
  virtual void Run();

  /// Set sigma
  SetMacro(Sigma, double);

  /// Get sigma
  GetMacro(Sigma, double);

};

#include <irtkGaussianBlurringWithPadding.h>

#endif
