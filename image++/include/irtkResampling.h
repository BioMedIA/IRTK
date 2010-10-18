/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKRESAMPLING_H

#define _IRTKRESAMPLING_H

#include <irtkImageToImage.h>

#include <irtkImageFunction.h>

#ifdef HAS_TBB

template <class VoxelType> class irtkMultiThreadedResampling;

#endif

/**
 * Class for resampling of images
 *
 * This class defines and implements the resampling of images with arbitrary
 * voxel dimensions.  The new image intensity of the voxels is calculated by
 * interpolation of the old image intensities. Possible interpolation schemes
 * are nearest neighbor, linear, cubic spline and B-spline interpolation.
 */

template <class VoxelType> class irtkResampling : public irtkImageToImage<VoxelType>
{

#ifdef HAS_TBB

  friend class irtkMultiThreadedResampling<VoxelType>;

#endif

protected:

  /// Voxel size of output after resampling
  double _XSize;
  double _YSize;
  double _ZSize;

  /// Interpolation used to interpolate output
  irtkImageFunction *_Interpolator;

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Initialize the filter
  virtual void Initialize();

public:

  /// Constructor
  irtkResampling(double, double, double);

  /// Set voxel size after resampling
  SetMacro(XSize, double);

  /// Get voxel size after resampling
  GetMacro(XSize, double);

  /// Set voxel size after resampling
  SetMacro(YSize, double);

  /// Get voxel size after resampling
  GetMacro(YSize, double);

  /// Set voxel size after resampling
  SetMacro(ZSize, double);

  /// Get voxel size after resampling
  GetMacro(ZSize, double);

  /// Set interpolator
  SetMacro(Interpolator, irtkImageFunction *);

  /// Run the resampling filter
  virtual void   Run();

};

#include <irtkResamplingWithPadding.h>

#endif
