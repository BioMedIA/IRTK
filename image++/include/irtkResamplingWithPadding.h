/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKRESAMPLINGWITHPADDING_H

#define _IRTKRESAMPLINGWITHPADDING_H

#ifdef HAS_TBB

template <class VoxelType> class irtkMultiThreadedResamplingWithPadding;

#endif

/**
 * Class for resampling of padded images
 *
 * This class defines and implements the resampling of images with arbitrary
 * voxel dimensions.  The new image intensity of the voxels is calculated by
 * interpolation of the old image intensities.  Only linear interpolation is
 * currently supported. If more than 50% of the voxels used for interpolation
 * have intensities smaller or equal to the padding value, the resampled
 * voxel will be filled with the padding value.
 */


template <class VoxelType> class irtkResamplingWithPadding : public irtkResampling<VoxelType>
{

#ifdef HAS_TBB

  friend class irtkMultiThreadedResamplingWithPadding<VoxelType>;

#endif

protected:

  /// Padding value
  VoxelType _PaddingValue;

  /// Returns whether the filter requires buffering
  virtual Bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Initialize the filter
  virtual void Initialize();

public:

  /// Constructor
  irtkResamplingWithPadding(double, double, double, VoxelType);

  /// Run the resampling filter
  virtual void Run();

  /// Set padding value
  SetMacro(PaddingValue, VoxelType);

  /// Get padding value
  GetMacro(PaddingValue, VoxelType);

};

#endif
