/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKNOISE_H

#define _IRTKNOISE_H

// irtk includes
#include <irtkImageToImage.h>

/**
 * Class for adding noise to images
 *
 * This class defines an abstract filter for adding noise to images.
 *
 */
template <class VoxelType> class irtkNoise : public irtkImageToImage<VoxelType>
{

protected:

  /// Initialization for the noise filter
  long int _Init;

  /// Amplitude for the noise filter
  double _Amplitude;

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor
  irtkNoise(double amplitude = 1);

  /// Destructor (empty).
  ~irtkNoise() {};

  /// Set amplitude
  SetMacro(Amplitude, double);

  /// Get amplitude
  GetMacro(Amplitude, double);

};

#include <irtkUniformNoise.h>
#include <irtkUniformNoiseWithPadding.h>
#include <irtkGaussianNoise.h>
#include <irtkGaussianNoiseWithPadding.h>
#include <irtkRicianNoise.h>
#include <irtkRicianNoiseWithPadding.h>

#endif
