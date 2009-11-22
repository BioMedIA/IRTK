/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKUNIFORMNOISEWITHPADDING_H

#define _IRTKUNIFORMNOISEWITHPADDING_H

#include <irtkUniformNoise.h>

/**
 * Class for adding uniform noise to images which are padded
 *
 * This class implements a filter for adding uniform noise to images which
 * are padded.
 *
 */

template <class VoxelType> class irtkUniformNoiseWithPadding : public irtkUniformNoise<VoxelType>
{

protected:
  /// Returns the name of the class.
  virtual const char* NameOfClass();

  /// The padding value.
  VoxelType _PaddingValue;

public:

  // Default constructor
  irtkUniformNoiseWithPadding();

  /** Constructor. Sets standard deviation of the noise distribution and
   *  the padding value.
   */
  irtkUniformNoiseWithPadding(double Amplitude, VoxelType PaddingValue);

  /// Destructor (empty).
  ~irtkUniformNoiseWithPadding() {};

  /// Set padding value
  SetMacro(PaddingValue, VoxelType);

  /// Get padding value
  GetMacro(PaddingValue, VoxelType);

  /// Run uniform noise filter
  virtual double Run(int, int, int, int);

};

#endif
