





#ifndef _IRTKGAUSSIANBLURRINGWITHPADDING2D_H

#define _IRTKGAUSSIANBLURRINGWITHPADDING2D_H

/**
 * Class for Gaussian blurring of padded images
 *
 * This class defines and implements the Gaussian blurring of padded images.
 * It takes 2D and 3D images but blurres only in the x and y direction.
 * The blurring is implemented by two successive 1D convolutions with a 1D
 * Gaussian kernel. If more than 50% of the voxels used for the convolution
 * have intensities smaller or equal to the padding value, the blurred voxel
 * will be filled with the padding value.
 */

template <class VoxelType> class irtkGaussianBlurringWithPadding2D : public irtkGaussianBlurring2D<VoxelType>
{

protected:

  /// Padding value
  VoxelType _PaddingValue;

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor
  irtkGaussianBlurringWithPadding2D(double, VoxelType);

  /// Run Gaussian blurring
  virtual void Run();

  /// Set padding value
  SetMacro(PaddingValue, VoxelType);

  /// Get padding value
  GetMacro(PaddingValue, VoxelType);

};

#endif
