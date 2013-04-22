/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : 
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : 
  Version   : 
  Changes   : $Author: bkainz $
  Original reduction implementation by: Noel Lopes, GPUMLib

=========================================================================*/


#pragma once

#include "cutil_math.h"
#include "irtkCUImage.h"

inline int divup(int a, int b) { return (a % b != 0) ? (a / b + 1) : (a / b); }
inline dim3 divup( uint2 a, dim3 b) { return dim3(divup(a.x, b.x), divup(a.y, b.y)); }
inline dim3 divup( dim3 a, dim3 b) { return dim3(divup(a.x, b.x), divup(a.y, b.y), divup(a.z, b.z)); }


// declare your kernel functions here!

inline __device__ float2 toFloat( const short2 & data ){
    return make_float2(data.x / 32766.0f, data.y);
}

inline __device__ float3 toFloat( const short3 & data ){
	return make_float3(data.x / 32766.0f, data.y, data.z);
}

inline __device__ short2 fromFloat( const float2 & data ){
    return make_short2(data.x * 32766.0f, data.y);
}

inline __device__ short3 fromFloat( const float3 & data ){
	return make_short3(data.x * 32766.0f, data.y, data.z);
}

//declarations
template <class T> 
void mask_image_gpu(int size, T *d_idata, T *d_odata, T *d_mdata, T nullVal, T maskValue);

//declarations
template <class T> 
void normalize_image_gpu(int size, T *d_idata, T *d_odata, T minVal, T maxVal);


#define INSTANCE_MACRO(TYPE) template void normalize_image_gpu<TYPE>(int size, TYPE *d_idata, TYPE *d_odata, TYPE minVal, TYPE maxVal);
INSTANCE_MACRO(char)
INSTANCE_MACRO(unsigned char)
INSTANCE_MACRO(short)
INSTANCE_MACRO(unsigned short)
INSTANCE_MACRO(int)
INSTANCE_MACRO(unsigned int)
INSTANCE_MACRO(float)
INSTANCE_MACRO(double)

#define INSTANCE_MACRO1(TYPE) template void mask_image_gpu<TYPE>(int size, TYPE *d_idata, TYPE *d_odata, TYPE *d_mdata, TYPE nullVal, TYPE maskValue);
INSTANCE_MACRO1(char)
INSTANCE_MACRO1(unsigned char)
INSTANCE_MACRO1(short)
INSTANCE_MACRO1(unsigned short)
INSTANCE_MACRO1(int)
INSTANCE_MACRO1(unsigned int)
INSTANCE_MACRO1(float)
INSTANCE_MACRO1(double)
