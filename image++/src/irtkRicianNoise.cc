/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

// irtk includes
#include <irtkImage.h>
#include <irtkNoise.h>

// boost includes
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

template <class VoxelType> irtkRicianNoise<VoxelType>::irtkRicianNoise(double Amplitude) : irtkNoise<VoxelType>(Amplitude)
{
}

template <class VoxelType> const char *irtkRicianNoise<VoxelType>::NameOfClass()
{
  return "irtkRicianNoise";
}

template <class VoxelType> double irtkRicianNoise<VoxelType>::Run(int x, int y, int z, int t)
{
  boost::mt19937 rng;
  rng.seed(this->_Init);

  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<boost::mt19937&,
			   boost::normal_distribution<> > var_nor(rng, nd);

#ifdef OLD_RICIAN
  // For backward consistency:
  if (this->_input->Get(x, y, z) <= 3*this->_Amplitude) {
    // Rayleigh distributed noise in background and air
    double s1 = var_nor();
    double s2 = var_nor();
    double M  = sqrt(s1*s1+s2*s2);
    return double(this->_input->Get(x, y, z, y)) + this->_Amplitude * M;
  } else {
    // Gaussian distributed noise in foreground
    return double(this->_input->Get(x, y, z, y)) + this->_Amplitude * var_nor();
  }

#else
  // Add Rician noise by treating the magnitude image as the real part (r) of a new complex variable (r,i),
  // adding Gaussian noise to the separate parts of this complex variable, before taking the magnitude again.
  // Thanks to Ged Ridgway <gerard.ridgway@ucl.ac.uk> for this contribution.
  double r = double(this->_input->Get(x, y, z, t)) + this->_Amplitude * var_nor();
  double i = this->_Amplitude * var_nor();
  return sqrt(r*r + i*i);

#endif
}

template class irtkRicianNoise<irtkBytePixel>;
template class irtkRicianNoise<irtkGreyPixel>;
template class irtkRicianNoise<irtkRealPixel>;

