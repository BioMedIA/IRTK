





#ifndef _IRTKGAUSSIANBLURRING2D_H

#define _IRTKGAUSSIANBLURRING2D_H

#include <irtkImageToImage.h>

/**
 * Class for Gaussian blurring of images
 *
 * This class defines and implements the Gaussian blurring of images. It takes
 * 2D and 3D images but blurres only in the x and y direction. The
 * blurring is implemented by two successive 1D convolutions with a 1D
 * Gaussian kernel.
 */

template <class VoxelType> class irtkGaussianBlurring2D : public irtkImageToImage<VoxelType>
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
  irtkGaussianBlurring2D(double);

  /// Destructor
  ~irtkGaussianBlurring2D();

  /// Run Gaussian blurring
  virtual void Run();

  /// Set sigma
  SetMacro(Sigma, double);

  /// Get sigma
  GetMacro(Sigma, double);

};

#include <irtkGaussianBlurringWithPadding2D.h>

#endif
