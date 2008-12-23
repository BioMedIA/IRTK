/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCONVOLUTION_H

#define _IRTKCONVOLUTION_H

#include <irtkImageToImage.h>

template <class VoxelType> class irtkConvolution : public irtkImageToImage<VoxelType>
{

protected:

  /// Flag whether to normalize convolution
  Bool _Normalization;

public:

  /// Constructor
  irtkConvolution(Bool = False);

  /// Set normalization on/off
  SetMacro(Normalization, Bool);

  /// Set normalization on/off
  GetMacro(Normalization, Bool);

};

// Convolution filters without padding
#include <irtkConvolution_1D.h>
#include <irtkConvolution_2D.h>
#include <irtkConvolution_3D.h>

// Convolution filters with padding
#include <irtkConvolutionWithPadding_1D.h>
#include <irtkConvolutionWithPadding_2D.h>
#include <irtkConvolutionWithPadding_3D.h>

#endif
