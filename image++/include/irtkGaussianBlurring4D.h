/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGAUSSIANBLURRING4D_H

#define _IRTKGAUSSIANBLURRING4D_H

#include <irtkImageToImage.h>

/**
 * Class for Gaussian blurring of image sequences
 *
 * This class defines and implements the Gaussian blurring of image sequences.
 * The blurring is implemented by four successive 1D convolutions with a 1D
 * Gaussian kernel.
 */

template <class VoxelType> class irtkGaussianBlurring4D : public irtkImageToImage<VoxelType>
{

protected:

  /// Sigma (standard deviation of Gaussian kernel)
  double _Sigma;

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor
  irtkGaussianBlurring4D(double);

  /// Destructor
  ~irtkGaussianBlurring4D();

  /// Run Gaussian blurring
  virtual void Run();

  /// Set sigma
  SetMacro(Sigma, double);

  /// Get sigma
  GetMacro(Sigma, double);

};

#endif

