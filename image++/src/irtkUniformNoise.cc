/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkNoise.h>

#include <nr.h>
#include <nrutil.h>

template <class VoxelType> irtkUniformNoise<VoxelType>::irtkUniformNoise(double Amplitude) : irtkNoise<VoxelType>(Amplitude)
{
  // Initialise ran2 with negative argument.
  long temp = -1 * this->_Init;
  (void) ran2(&temp);
}

template <class VoxelType> const char *irtkUniformNoise<VoxelType>::NameOfClass()
{
  return "irtkUniformNoise";
}

template <class VoxelType> double irtkUniformNoise<VoxelType>::Run(int x, int y, int z, int t)
{
  return double(this->_input->Get(x, y, z, t)) + this->_Amplitude * ran2(&this->_Init);
}

template class irtkUniformNoise<irtkBytePixel>;
template class irtkUniformNoise<irtkGreyPixel>;
template class irtkUniformNoise<irtkRealPixel>;
