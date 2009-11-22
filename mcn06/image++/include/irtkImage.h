/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGE_H

#define _IRTKIMAGE_H

// Basic voxel classes
#include <irtkVoxel.h>

// Basic image class
template <class Type> class irtkGenericImage;

// Includes
#include <irtkGeometry.h>

#include <irtkBaseImage.h>
#include <irtkGenericImage.h>

/// Unsigned char image
typedef class irtkGenericImage<irtkBytePixel> irtkByteImage;
/// Short image
typedef class irtkGenericImage<irtkGreyPixel> irtkGreyImage;
/// Float image
typedef class irtkGenericImage<irtkRealPixel> irtkRealImage;

#ifndef _IMPLEMENTS_GENERICIMAGE_

extern template class irtkGenericImage<irtkBytePixel>;
extern template class irtkGenericImage<irtkGreyPixel>;
extern template class irtkGenericImage<irtkRealPixel>;

#endif

#endif
