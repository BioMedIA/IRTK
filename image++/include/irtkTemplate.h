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

#ifdef _IMPLEMENTS_GENERICIMAGE_

template class irtkGenericImage<char>;
template class irtkGenericImage<unsigned char>;
template class irtkGenericImage<short>;
template class irtkGenericImage<unsigned short>;
template class irtkGenericImage<int>;
template class irtkGenericImage<unsigned int>;
template class irtkGenericImage<float>;
template class irtkGenericImage<double>;

template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<float> &);

template irtkGenericImage<float>& irtkGenericImage<float>::operator=(const irtkGenericImage<short> &);
template irtkGenericImage<float>& irtkGenericImage<float>::operator=(const irtkGenericImage<double> &);
template irtkGenericImage<short>& irtkGenericImage<short>::operator=(const irtkGenericImage<float> &);
template irtkGenericImage<short>& irtkGenericImage<short>::operator=(const irtkGenericImage<double> &);
template irtkGenericImage<double>& irtkGenericImage<double>::operator=(const irtkGenericImage<float> &);
template irtkGenericImage<double>& irtkGenericImage<double>::operator=(const irtkGenericImage<short> &);

#endif
