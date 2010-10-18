/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCONVOLUTIONWITHPADDING_1D_H

#define _IRTKCONVOLUTIONWITHPADDING_1D_H

/**
 * Class for one-dimensional convolution of padded images.
 *
 * This class defines and implements one-dimensional convolutions of a padded
 * image with a filter kernel. The convolution is computed along the x-axis.
 * This class assumes that the filter kernel is one-dimensional and its size
 * along the y- and z-axis must be 1. The filter ignores any padded values in
 * the image.
 */

template <class VoxelType> class irtkConvolutionWithPadding_1D : public irtkConvolution_1D<VoxelType>
{

protected:

  /** Padding value. The padding value defines which voxels are ignored
   *  during the calculation of the convolution. All voxels whose value
   *  is smaller or equal to the padding value are ignored.
   */
  VoxelType _padding;

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /** Run the convolution filter. This method is protected and should only
   *  be called from within public member function Run(). This method also
   *  overrides the member function Run() of the base class and ignores any
   *  padded values during the convolution. If more than 50\% of the voxels
   *  are padded, it returns the padding value, otherwise it returns the
   *  result of the convolution.
   */
  virtual double Run(int, int, int, int);

public:

  /// Constructor
  irtkConvolutionWithPadding_1D(VoxelType, bool = false);

  /// Put the padding value
  void      PutPaddingValue(VoxelType);

  /// Get the padding value
  VoxelType GetPaddingValue();
};

#endif
