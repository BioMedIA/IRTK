/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    :
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      :
  Version   :
  Changes   : $Author: bkainz $

=========================================================================*/

#ifndef NULLPTR_H_
#define NULLPTR_H_


// keywords required if nullptr used in CUDA code (compiled by nvcc)
#ifdef __CUDACC__
#  include <cuda.h>
#  define _NULLPTR_DECL __host__ __device__
#else
#  define _NULLPTR_DECL
#endif


const                                      // this is a const object...
const class {
public:
   template<class T>                       // convertible to any type
   _NULLPTR_DECL operator T*() const       // of null non-member
    { return 0; }                          // pointer...
   template<class C, class T>              // or any type of null
   _NULLPTR_DECL operator T C::*() const   // member pointer...
    { return 0; }
private:
  _NULLPTR_DECL void operator&() const;    // whose address can't be taken
} nullptr = {};


#undef _NULLPTR_DECL

#endif // NULLPTR_H_
