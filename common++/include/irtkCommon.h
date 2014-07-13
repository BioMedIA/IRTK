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

#ifndef _IRTKCOMMON_H

#define _IRTKCOMMON_H

// C++ header files
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <algorithm>
#include <string>
#include <limits>

// C header files
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>

#ifdef HAS_ZLIB
#include <zlib.h>
#endif

#ifdef WIN32

// Windows called copysign differently before C++11
#define copysign _copysign

// Windows is missing M_PI constants
#define M_PI		3.14159265358979323846

#define NOMINMAX

// Windows specific header file
#include <windows.h>

#endif

// Use standard namespace
using namespace std;

#define SetMacro(name,type) \
void Set##name (type arg) { this->_##name = arg; }

#define GetMacro(name,type) \
type Get##name () { return this->_##name; }

extern const char *dirname2(const char *path);
extern char       *basename2(const char *filename);

extern void PrintVersion(ostream &, const char*);

extern int   ReadInt   (ifstream &in);
extern float ReadFloat (ifstream &in);
extern char *ReadString(ifstream &in);

#ifdef HAS_ZLIB
extern int   ReadCompressed(gzFile, char *, long, long);
#else
extern int   ReadCompressed(FILE *, char *, long, long);
#endif

#define round round2

inline int round(double x)
{
  return x > 0 ? int(x + 0.5) : int(x - 0.5);
}

extern void swap16(char *, char *, long);
extern void swap32(char *, char *, long);
extern void swap64(char *, char *, long);

class sort_indices
{
  private:
    float* _array;

  public:
    sort_indices(float* array) : _array(array) {}
    bool operator()(int i, int j) { return _array[i] < _array[j]; } 
};

extern double weightedmedian(int, double,double, float*, float*);

// Orientation codes (same as NIFTI)
#define IRTK_L2R  1    /* Left to Right         */
#define IRTK_R2L  2    /* Right to Left         */
#define IRTK_P2A  3    /* Posterior to Anterior */
#define IRTK_A2P  4    /* Anterior to Posterior */
#define IRTK_I2S  5    /* Inferior to Superior  */
#define IRTK_S2I  6    /* Superior to Inferior  */

#include <irtkObject.h>
#include <irtkCifstream.h>
#include <irtkCofstream.h>
#include <irtkAllocate.h>
#include <irtkDeallocate.h>
#include <irtkException.h>
#include <irtkParallel.h>

#ifdef HAS_VTK
#  include <vtkConfigure.h>
#endif


#endif
