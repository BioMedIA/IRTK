/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGAUSSIANNOISEWITHPADDING_H

#define _IRTKGAUSSIANNOISEWITHPADDING_H

#include <irtkGaussianNoise.h>

/**
 * Class for adding Gaussian noise to images which are padded
 *
 * This class implements a filter for adding Gaussian noise to images which
 * are padded.
 *
 */

template <class VoxelType> class irtkGaussianNoiseWithPadding : public irtkGaussianNoise<VoxelType>
{

protected:
  /// Returns the name of the class.
  virtual const char* NameOfClass();

  /// The padding value.
  VoxelType _PaddingValue;

public:

  // Default constructor
  irtkGaussianNoiseWithPadding();

  /** Constructor.  Sets mean value and standard deviation of the noise
   *  distribution, the padding value, and the image range.
   */
  irtkGaussianNoiseWithPadding(double Mean, double Sigma, VoxelType PaddingValue, VoxelType MinVal, VoxelType MaxVal);

  /// Destructor (empty).
  ~irtkGaussianNoiseWithPadding() {};

  /// Set padding value
  SetMacro(PaddingValue, VoxelType);

  /// Get padding value
  GetMacro(PaddingValue, VoxelType);

  /// Adds noise to input image.
  virtual double Run(int, int, int, int);

};

#endif
