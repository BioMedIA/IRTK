/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKRICIANNOISE_H

#define _IRTKRICIANNOISE_H

/**
 * Class for adding Rician noise to images
 *
 * This class implements a filter for adding Rician noise to images. The
 * Rician distribution is approiximately Gaussian for high intensity
 * signals, and Rayleigh distributed for low intensities. The Rayleigh
 * intensity distribution can be expressed as:
 * $P_M(M)=\frac{m}{\sigma^2}exp\left(\frac{-M^2}{2\sigma^2}\right)
 * where $M$ is the actual intensity and $\sigma$ is the standard
 * deviation of Gaussian noise. For more information, see Holden et al.,
 * IEEE-TMI 19(2) 2000.
 *
 */

template <class VoxelType> class irtkRicianNoise : public irtkNoise<VoxelType>
{

protected:

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor
  irtkRicianNoise(double amplitude = 1);

  /// Destructor (empty).
  ~irtkRicianNoise() {};

  /// Run rician noise filter
  virtual double Run(int, int, int, int);

};

#endif

