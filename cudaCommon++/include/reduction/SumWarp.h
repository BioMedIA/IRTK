/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : 
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : 
  Version   : 
  Changes   : $Author: bkainz $
  GPL licensed file: Original reduction implementation by: Noel Lopes, GPUMLib

=========================================================================*/


#ifndef GPUMLib_SumWarp_h
#define GPUMLib_SumWarp_h

#include "../reduction/CudaDefinitions.h"

template <class VoxelType, int blockSize> __device__ __forceinline__ void SumBeforeWarp(VoxelType * s) {
	if (blockSize >= 1024) {
		if (threadIdx.x < 512) s[threadIdx.x] += s[threadIdx.x + 512];
		__syncthreads();
	}

	if (blockSize >= 512) {
		if (threadIdx.x < 256) s[threadIdx.x] += s[threadIdx.x + 256];
		__syncthreads();
	}

	if (blockSize >= 256) {
		if (threadIdx.x < 128) s[threadIdx.x] += s[threadIdx.x + 128];
		__syncthreads();
	}

	if (blockSize >= 128) {
		if (threadIdx.x < 64) s[threadIdx.x] += s[threadIdx.x + 64];
		__syncthreads();
	}
}

template <class VoxelType, int blockSize> __device__ __forceinline__ void SumWarp(volatile VoxelType * s) {
	if (blockSize >= 64) s[threadIdx.x] += s[threadIdx.x + 32];
	if (blockSize >= 32) s[threadIdx.x] += s[threadIdx.x + 16];
	if (blockSize >= 16) s[threadIdx.x] += s[threadIdx.x + 8];
	if (blockSize >= 8) s[threadIdx.x] += s[threadIdx.x + 4];
	if (blockSize >= 4) s[threadIdx.x] += s[threadIdx.x + 2];
	if (blockSize >= 2) s[threadIdx.x] += s[threadIdx.x + 1];
}

#endif