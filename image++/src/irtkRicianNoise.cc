/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

// irtk includes
#include <irtkImage.h>
#include <irtkNoise.h>

// recipe includes
#include <nr.h>
#include <nrutil.h>

template <class VoxelType> irtkRicianNoise<VoxelType>::irtkRicianNoise(double Amplitude) : irtkNoise<VoxelType>(Amplitude)
{
}

template <class VoxelType> const char *irtkRicianNoise<VoxelType>::NameOfClass()
{
  return "irtkRicianNoise";
}

template <class VoxelType> double irtkRicianNoise<VoxelType>::Run(int x, int y, int z, int t)
{
#ifdef OLD_RICIAN
  // For bug-ward consistency:

  if (this->_input->Get(x, y, z) <= 3*this->_Amplitude) {
    // Rayleigh distributed noise in background and air
    double s1 = gasdev(&this->_Init);
    double s2 = gasdev(&this->_Init);
    double M  = sqrt(s1*s1+s2*s2);
    return double(this->_input->Get(x, y, z, y)) + this->_Amplitude * M;
  } else {
    // Gaussian distributed noise in foreground
    return double(this->_input->Get(x, y, z, y)) + this->_Amplitude * gasdev(&this->_Init);
  }

#else
  // Add Rician noise by treating the magnitude image as the real part (r) of a new complex variable (r,i),
  // adding Gaussian noise to the separate parts of this complex variable, before taking the magnitude again.
  // Thanks to Ged Ridgway <gerard.ridgway@ucl.ac.uk> for this contribution.
  double r = double(this->_input->Get(x, y, z, t)) + this->_Amplitude * gasdev(&this->_Init);
  double i = this->_Amplitude * gasdev(&this->_Init);
  return sqrt(r*r + i*i);

#endif
}

template class irtkRicianNoise<irtkBytePixel>;
template class irtkRicianNoise<irtkGreyPixel>;
template class irtkRicianNoise<irtkRealPixel>;

