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


#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <math.h>
#include <limits>
#include <vector_types.h>
#include <vector_functions.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <math_functions.h>
#include <math_constants.h>

#include "irtkCUImage_kernels.cuh"

#if (__CUDA_ARCH__ < 200)
#define int_mult(x,y)	__mul24(x,y)	
#else
#define int_mult(x,y)	x*y
#endif


//TODO add some image statistics functions here if necessary

// simple filter only evaluated array w/o class overhead me transfers
template <class T>
__global__ void
normalize_kernel(T *g_idata, T *g_odata, T minv, T maxv, T typemin, T typemax)
{
	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
	g_odata[i] = (T)(typemin + (typemax-typemin)*((double)g_idata[i] - (double)minv) /((double)maxv - (double)minv));
}

//specialization for float and double -> normalization from 0 to 1
template <>
__global__ void
normalize_kernel<float>(float *g_idata, float *g_odata, float minv, float maxv, float typemin, float typemax)
{
	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
	g_odata[i] = (float)(((double)g_idata[i] - (double)minv) /((double)maxv - (double)minv));
}

template <>
__global__ void
normalize_kernel<double>(double *g_idata, double *g_odata, double minv, double maxv, double typemin, double typemax)
{
	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
	g_odata[i] = (((double)g_idata[i] - (double)minv) /((double)maxv - (double)minv));
}


// simple filter only evaluated array w/o class overhead me transfers
template <class T>
__global__ void
mask_image_kernel(T *g_idata, T *g_odata, T *g_mdata, T nullVal, T maskValue)
{
	unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(g_mdata[i] == maskValue)
		g_odata[i] = nullVal;
	else
		g_odata[i] = g_idata[i];
}


template <class T>
void
mask_image_gpu(int size, T *d_idata, T *d_odata, T *d_mdata, T nullVal, T maskValue )
{
	int threadsPerBlock = 256;

    dim3 dimBlock(threadsPerBlock, 1, 1);
    dim3 dimGrid((size + threadsPerBlock - 1) / threadsPerBlock, 1, 1);

	mask_image_kernel<T><<<dimGrid, dimBlock>>>(d_idata, d_odata, d_mdata, nullVal, maskValue);
	cudaDeviceSynchronize();
}

template <class T>
void
normalize_image_gpu(int size, T *d_idata, T *d_odata, T minVal, T maxVal)
{
	int threadsPerBlock = 512;

    dim3 dimBlock(threadsPerBlock, 1, 1);
    dim3 dimGrid((size + threadsPerBlock - 1) / threadsPerBlock, 1, 1);

	normalize_kernel<T><<<dimGrid, dimBlock>>>(d_idata, d_odata, minVal, maxVal, numeric_limits<T>::min(), numeric_limits<T>::max());
	cudaDeviceSynchronize();
}


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
