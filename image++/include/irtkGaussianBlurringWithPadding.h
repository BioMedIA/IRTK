/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGAUSSIANBLURRINGWITHPADDING_H

#define _IRTKGAUSSIANBLURRINGWITHPADDING_H

/**
 * Class for Gaussian blurring of padded images
 *
 * This class defines and implements the Gaussian blurring of padded images.
 * The blurring is implemented by three successive 1D convolutions with a 1D
 * Gaussian kernel. If more than 50% of the voxels used for the convolution
 * have intensities smaller or equal to the padding value, the blurred voxel
 * will be filled with the padding value.
 */

template <class VoxelType> class irtkGaussianBlurringWithPadding : public irtkGaussianBlurring<VoxelType>
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
  irtkGaussianBlurringWithPadding(double, VoxelType);

  /// Run Gaussian blurring
  virtual void Run();

  /// Set padding value
  SetMacro(PaddingValue, VoxelType);

  /// Get padding value
  GetMacro(PaddingValue, VoxelType);

};

#endif
