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

#include <irtkConvolution.h>

template <class VoxelType> irtkConvolution<VoxelType>::irtkConvolution(Bool Normalization)
{
  _Normalization = Normalization;
}

template class irtkConvolution<unsigned char>;
template class irtkConvolution<short>;
template class irtkConvolution<unsigned short>;
template class irtkConvolution<float>;
template class irtkConvolution<double>;
