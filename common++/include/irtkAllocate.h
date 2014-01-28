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

#ifndef _IRTKALLOCATE_H

#define _IRTKALLOCATE_H

#include <irtkCommon.h>

/// Allocate 1-dimensional array and initialize elements
///
/// @code
/// // Allocate an array of 10 pointers to a grey image and initialize them to NULL
/// irtkGreyImage *images = Allocate<irtkGreyImage *>(10);
/// // Allocate an array of 5 integer values and initialize them to -1
/// int *values = Allocate<int>(5, -1);
/// @endcode
template <typename Type> inline Type *Allocate(int n, const Type &init = Type())
{
  Type *p = new Type[n];
  for (int i = 0; i < n; i++) p[i] = init;
  return p;
}

template <class Type> inline Type **Allocate(Type **matrix, int x, int y)
{
  int i;

  if ((matrix = new Type *[y]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << "\n";
    exit(1);
  }

  if ((matrix[0] = new Type[y*x]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << "\n";
    exit(1);
  }

  for (i = 1; i < y; i++) {
    matrix[i] = matrix[i-1] + x;
  }

#ifdef DEBUG
  memory_allocated += x*y*sizeof(Type);
  cout << "Allocate: Memory allocated is " << memory_allocated << " bytes.\n";
#endif

  return matrix;
}

template <class Type> inline Type ***Allocate(Type ***matrix, int x, int y, int z)
{
  int i, j;

  if ((matrix = new Type **[z]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
    cerr << z << "\n";
    exit(1);
  }

  if ((matrix[0] = new Type *[z*y]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
    cerr << z << "\n";
    exit(1);
  }

  for (i = 1; i < z; i++) {
    matrix[i] = matrix[i-1] + y;
  }

  if ((matrix[0][0] = new Type[z*y*x]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
    cerr << z << "\n";
    exit(1);
  }

  for (i = 0; i < z; i++) {
    for (j = 0; j < y; j++) {
      matrix[i][j] = matrix[0][0] + (i*y+j)*x;
    }
  }

#ifdef DEBUG
  memory_allocated += x*y*z*sizeof(Type);
  cout << "Allocate: Memory allocated is " << memory_allocated << " bytes.\n";
#endif

  return matrix;
}

template <class Type> inline Type ****Allocate(Type ****matrix, int x, int y, int z, int t)
{
  int i, j, k;

  if ((matrix = new Type ***[t]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
    cerr << z << " x " << t << "\n";
    exit(1);
  }

  if ((matrix[0] = new Type **[t*z]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
    cerr << z << " x " << t << "\n";
    exit(1);
  }

  for (i = 1; i < t; i++) {
    matrix[i] = matrix[i-1] + z;
  }

  if ((matrix[0][0] = new Type*[t*z*y]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
    cerr << z << " x " << t << "\n";
    exit(1);
  }

  for (i = 0; i < t; i++) {
    for (j = 0; j < z; j++) {
      matrix[i][j] = matrix[0][0] + i*z*y + j*y;
    }
  }

  if ((matrix[0][0][0] = new Type[t*z*y*x]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x ";
    cerr << z << " x " << t << "\n";
    exit(1);
  }

  for (i = 0; i < t; i++) {
    for (j = 0; j < z; j++) {
      for (k = 0; k < y; k++) {
        matrix[i][j][k] = matrix[0][0][0] + i*z*y*x + j*y*x + k*x;
      }
    }
  }

#ifdef DEBUG
  memory_allocated += x*y*z*t*sizeof(Type);
  cout << "Allocate: Memory allocated is " << memory_allocated << " bytes.\n";
#endif

  return matrix;
}

#endif
