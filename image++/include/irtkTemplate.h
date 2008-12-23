/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifdef _IMPLEMENTS_GENERICIMAGE_

template class irtkGenericImage<irtkBytePixel>;
template class irtkGenericImage<irtkGreyPixel>;
template class irtkGenericImage<irtkRealPixel>;
template class irtkGenericImage<irtkRGBPixel>;
template class irtkGenericImage<irtkVector3D<char> >;
template class irtkGenericImage<irtkVector3D<short> >;
template class irtkGenericImage<irtkVector3D<float> >;
template class irtkGenericImage<irtkVector3D<double> >;

#endif
