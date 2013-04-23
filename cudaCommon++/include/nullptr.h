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

#include <cuda.h>

#ifdef __CUDACC__
const                        // this is a const object...
class {
public:
   template<class T>          // convertible to any type
   __host__ __device__ operator T*() const      // of null non-member
    { return 0; }            // pointer...
   template<class C, class T> // or any type of null
    __host__ __device__ operator T C::*() const  // member pointer...
    { return 0; }
private:
  __host__ __device__ void operator&() const;    // whose address can't be taken
} nullptr = {};

#endif 


#endif // NULLPTR_H_
