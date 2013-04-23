/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : 
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : 
  Version   : 
  Changes   : $Author: bkainz $

=========================================================================*/

#ifndef _IRTKCUIMAGE_H
#define _IRTKCUIMAGE_H

// Includes
#include <irtkVoxel.h>

// Includes
#include <irtkGeometry.h>

#include <irtkBaseImage.h>

#include <irtkCUGenericImage.h>
#include "irtkCUSharedLibMacros.h"

// Basic image class
template <class Type> class irtkCULib_DLLAPI irtkCUGenericImage;

/// Unsigned char image
typedef class irtkCULib_DLLAPI irtkCUGenericImage<irtkBytePixel> irtkCUByteImage;
/// Short image
typedef class irtkCULib_DLLAPI irtkCUGenericImage<irtkGreyPixel> irtkCUGreyImage;
/// Float image
typedef class irtkCULib_DLLAPI irtkCUGenericImage<irtkRealPixel> irtkCURealImage;

#ifndef _IMPLEMENTS_CUGENERICIMAGE_

extern template class irtkCULib_DLLAPI irtkCUGenericImage<irtkBytePixel>;
extern template class irtkCULib_DLLAPI irtkCUGenericImage<irtkGreyPixel>;
extern template class irtkCULib_DLLAPI irtkCUGenericImage<irtkRealPixel>;

#endif

#ifdef HAS_OPENCV

#include <irtkImageToOpenCv.h>

#endif

#endif
