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

#include <irtkImage.h>

#include <irtkNoise.h>

#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>

template <class VoxelType> irtkUniformNoise<VoxelType>::irtkUniformNoise(double Amplitude) : irtkNoise<VoxelType>(Amplitude)
{
  boost::mt19937 rng;
  rng.seed(this->_Init);

  boost::uniform_int<> ud(0, 1);
  boost::variate_generator<boost::mt19937&,
			   boost::uniform_int<> > uni_engine(rng, ud); 
  (void) uni_engine();
}

template <class VoxelType> const char *irtkUniformNoise<VoxelType>::NameOfClass()
{
  return "irtkUniformNoise";
}

template <class VoxelType> double irtkUniformNoise<VoxelType>::Run(int x, int y, int z, int t)
{
  boost::mt19937 rng;
  rng.seed(this->_Init);

  boost::uniform_int<> ud(0, this->_Amplitude);
  boost::variate_generator<boost::mt19937&,
			   boost::uniform_int<> > uni_engine(rng, ud);
  return double(this->_input->Get(x, y, z, t)) + uni_engine();
}

template class irtkUniformNoise<irtkBytePixel>;
template class irtkUniformNoise<irtkGreyPixel>;
template class irtkUniformNoise<irtkRealPixel>;
