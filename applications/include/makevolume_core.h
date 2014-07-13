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

#ifndef makevolume_core_h
#define makevolume_core_h

// functor for indexing
template<class T>
struct index_sorting {
  const T _values;

public:
  index_sorting(const T values) : _values(values) {};
  bool operator() (long unsigned a, long unsigned b) const {
    return _values[a] < _values[b]; }    
};

long unsigned* indexing(long unsigned n, float* array);

#endif // makevolume_core_h
