/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <sys/types.h>

#ifndef WIN32
#include <sys/time.h>
#endif

#include <irtkImage.h>

#include <irtkNoise.h>

template <class VoxelType> irtkNoise<VoxelType>::irtkNoise(double Amplitude)
{
#ifndef WIN32
  timeval tv;
  gettimeofday(&tv, NULL);
  _Init = tv.tv_usec;
  _Amplitude = Amplitude;
#else
  cerr << "irtkNoise<VoxelType>::irtkNoise: Not implemented yet for Windows" << endl;
  exit(1);
#endif
}

template <class VoxelType> Bool irtkNoise<VoxelType>::RequiresBuffering()
{
  return False;
}

template <class VoxelType> const char *irtkNoise<VoxelType>::NameOfClass()
{
  return "irtkNoise";
}

template class irtkNoise<irtkBytePixel>;
template class irtkNoise<irtkGreyPixel>;
template class irtkNoise<irtkRealPixel>;
