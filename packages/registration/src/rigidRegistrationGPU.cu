/*
 * Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation and 
 * any modifications thereto.  Any use, reproduction, disclosure, or distribution 
 * of this software and related documentation without an express license 
 * agreement from NVIDIA Corporation is strictly prohibited.
 * 
 */

/* Example of integrating CUDA functions into an existing 
 * application / framework.
 * Host part of the device code.
 * Compiled with Cuda compiler.
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <ctime>
// includes, kernels
#include <sm_11_atomic_functions.h>
#include "C:\ProgramData\NVIDIA Corporation\NVIDIA GPU Computing SDK\C\common\inc\cutil_inline.h"

#define SHORT_SIZE = 32767;
typedef unsigned int  uint;
typedef unsigned char uchar;

__shared__ float shared[32];
__shared__ int sharedSamples[256];

texture<short, 3, cudaReadModeElementType> target;  // 3D texture
texture<short, 3, cudaReadModeNormalizedFloat> source;  // 3D texture
cudaArray *d_sourceArray = 0; // volume array used
cudaArray *d_targetArray = 0; // volume array used
////////////////////////////////////////////////////////////////////////////////
// declaration, forward

////////////////////////////////////////////////////////////////////////////////
//! Entry point for Cuda functionality on host side
//! @param argc  command line argument count
//! @param argv  command line arguments
//! @param data  data to process on the device
//! @param len   len of \a data
////////////////////////////////////////////////////////////////////////////////
__global__ void
kernel_outputNormSourceVals(float *d_output, float *d_matrix,int imageX, int imageY, int imageZ, int modX, int modY)
{
		//Thread id within block
	uint tx = threadIdx.x;
    uint ty = threadIdx.y;

	uint bx = blockIdx.x;	
	uint by = blockIdx.y;	

	uint i = (bx % (gridDim.x / imageZ)) * blockDim.x + tx;
	uint j =  by * blockDim.y + ty;
	uint z = bx / (gridDim.x / imageZ); // ranges from 1 to 149

	//Calculate warpID

	shared[tx] = d_matrix[tx];
    int int_final_res = -1;

	__syncthreads();

	int target_value = tex3D(target, i, j, z);
	if ((target_value >= 0) && (i < imageX) && (j < imageY) && (z < imageZ)) {
        		  
          // Perform the matrix transformation on the current pixel
		  float _x = (shared[0+0] * i + shared[0+1] * j + shared[0+2] * z + shared[0+3]);
		  float _y = (shared[4+0] * i + shared[4+1] * j + shared[4+2] * z + shared[4+3]);
		  float _z = (shared[8+0] * i + shared[8+1] * j + shared[8+2] * z + shared[8+3]);

        		// Check whether transformed point is inside source volume
		if ((_x > 0) && (_x < imageX-1) &&
			(_y > 0) && (_y < imageY-1) &&
			(_z > 0) && (_z < imageZ-1)) {
            
			float source_value = tex3D(source, i, j, z);
            d_output[i + j * imageX + z * imageX * imageY] = source_value * 32767;
        }
    }
}

__global__ void
kernel_atomicRegister(int *d_output,float *d_matrix, int imageX,int imageY, int imageZ,int modX, int modY)
{

	//Thread id within block, i.e. ranges from 1-16
	uint tx = threadIdx.x;
    uint ty = threadIdx.y;

	uint bx = blockIdx.x;	//ranges from 1 to 1043
	uint by = blockIdx.y;	//ranges from 1 to 20

	uint i = (bx % (gridDim.x / imageZ)) * blockDim.x + tx;
	uint j =  by * blockDim.y + ty;
	uint z = bx / (gridDim.x / imageZ); // ranges from 1 to 149

	shared[tx] = d_matrix[tx];
    int int_final_res = -1;

	__syncthreads();

	if ((tex3D(target, i, j, z) >= 0) && (i < imageX) && (j < imageY) && (z < imageZ)) {
        		  
          // Perform the matrix transformation on the current pixel
		  float _x = (shared[0+0] * i + shared[0+1] * j + shared[0+2] * z + shared[0+3]);
		  float _y = (shared[4+0] * i + shared[4+1] * j + shared[4+2] * z + shared[4+3]);
		  float _z = (shared[8+0] * i + shared[8+1] * j + shared[8+2] * z + shared[8+3]);

        		// Check whether transformed point is inside source volume
		if ((_x > 0) && (_x < imageX-1) &&
			(_y > 0) && (_y < imageY-1) &&
			(_z > 0) && (_z < imageZ-1)) {
            
			float svoxel = tex3D(source, _x, _y, _z);
            int_final_res = svoxel * 32767 * 64 + (tex3D(target, i, j, z));
			atomicAdd(d_output+(int_final_res),1);
        }
    }

}

__global__ void
kernel_registerSamples(int *d_output,int *d_samplesArray, float *d_matrix, int imageX,int imageY, int imageZ,int modX, int modY)
{
	//Thread id within block
	uint tx = threadIdx.x;
    uint ty = threadIdx.y;

	uint bx = blockIdx.x;	
	uint by = blockIdx.y;	

	uint i = (bx % (gridDim.x / imageZ)) * blockDim.x + tx;
	uint j =  by * blockDim.y + ty;
	uint z = bx / (gridDim.x / imageZ); // ranges from 1 to 149

	//Calculate warpID

	shared[tx] = d_matrix[tx];
    int int_final_res = -1;

	__syncthreads();

	if ((tex3D(target, i, j, z) >= 0) && (i < imageX) && (j < imageY) && (z < imageZ)) {
        		  
          // Perform the matrix transformation on the current pixel
		  float _x = (shared[0+0] * i + shared[0+1] * j + shared[0+2] * z + shared[0+3]);
		  float _y = (shared[4+0] * i + shared[4+1] * j + shared[4+2] * z + shared[4+3]);
		  float _z = (shared[8+0] * i + shared[8+1] * j + shared[8+2] * z + shared[8+3]);

        		// Check whether transformed point is inside source volume
		if ((_x > 0) && (_x < imageX-1) &&
			(_y > 0) && (_y < imageY-1) &&
			(_z > 0) && (_z < imageZ-1)) {
            
			float svoxel = tex3D(source, _x, _y, _z);
            int_final_res = svoxel * 32767 * 64 + (tex3D(target, i, j, z));
        }
    }


	if(int_final_res > -1){
		sharedSamples[tx + blockDim.x * ty] = 1;
	} else {
		sharedSamples[tx + blockDim.x * ty] = 0;
	}

	d_output[i + j * modX + z * modX * modY] = int_final_res;
	 __syncthreads();

	//Perform Reduction
    // do reduction in shared mem
    for(unsigned int s=256/2; s>0; s>>=1)  
    {
        if (tx + ty*blockDim.x < s) 
        {
            sharedSamples[tx + ty*blockDim.x] += sharedSamples[tx + ty*blockDim.x + s];
        }
        __syncthreads();
    }


	if(tx + ty * blockDim.x == 0) 
		d_samplesArray[blockIdx.x +  blockIdx.y * gridDim.x] =  sharedSamples[0];//sharedSamples[0];

}

__global__ void
kernel_register(int *d_output,float *d_matrix, int imageX,int imageY, int imageZ,int modX, int modY)
{
	//Thread id within block
	uint tx = threadIdx.x;
    uint ty = threadIdx.y;

	uint bx = blockIdx.x;	
	uint by = blockIdx.y;	

	uint i = (bx % (gridDim.x / imageZ)) * blockDim.x + tx;
	uint j =  by * blockDim.y + ty;
	uint z = bx / (gridDim.x / imageZ); // ranges from 1 to 149

	//Calculate warpID

	shared[tx] = d_matrix[tx];
    int int_final_res = -1;

	__syncthreads();

	int target_value = tex3D(target, i, j, z);
	if ((target_value >= 0) && (i < imageX) && (j < imageY) && (z < imageZ)) {
        		  
          // Perform the matrix transformation on the current pixel
		  float _x = (shared[0+0] * i + shared[0+1] * j + shared[0+2] * z + shared[0+3]);
		  float _y = (shared[4+0] * i + shared[4+1] * j + shared[4+2] * z + shared[4+3]);
		  float _z = (shared[8+0] * i + shared[8+1] * j + shared[8+2] * z + shared[8+3]);

        		// Check whether transformed point is inside source volume
		if ((_x > 0) && (_x < imageX-1) &&
			(_y > 0) && (_y < imageY-1) &&
			(_z > 0) && (_z < imageZ-1)) {
            
			float source_value = tex3D(source, _x, _y, _z) * 32767;
            int_final_res = source_value * 64 + target_value;
        }
    }

	d_output[i + j * modX + z * modX * modY] = int_final_res;
}

extern "C"
void initCuda(short *h_target, short *h_source, cudaExtent volumeSize)
{
    // create 3D arrays
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<short>();
    cutilSafeCall( cudaMalloc3DArray(&d_targetArray, &channelDesc, volumeSize) );
    cutilSafeCall( cudaMalloc3DArray(&d_sourceArray, &channelDesc, volumeSize) );

    // copy data to 3D array
    cudaMemcpy3DParms copyParams = {0};

	//TARGET
    copyParams.srcPtr   = make_cudaPitchedPtr((void*)h_target, volumeSize.width*sizeof(short), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_targetArray;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    cutilSafeCall( cudaMemcpy3D(&copyParams) );

	//SOURCE
    copyParams.srcPtr   = make_cudaPitchedPtr((void*)h_source, volumeSize.width*sizeof(short), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_sourceArray;
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;
    cutilSafeCall( cudaMemcpy3D(&copyParams) );

    // set texture parameters
    target.normalized = false;                      // access with normalized texture coordinates
    target.filterMode = cudaFilterModePoint;      // linear interpolation
    target.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates
    target.addressMode[1] = cudaAddressModeClamp;
    target.addressMode[2] = cudaAddressModeClamp;

    // bind array to 3D texture
    cutilSafeCall(cudaBindTextureToArray(target, d_targetArray, channelDesc));

    // set texture parameters
    source.normalized = false;                      // access with normalized texture coordinates
    source.filterMode = cudaFilterModePoint;      // linear interpolation
    source.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates
    source.addressMode[1] = cudaAddressModeClamp;
    source.addressMode[2] = cudaAddressModeClamp;

    // bind array to 3D texture
    cutilSafeCall(cudaBindTextureToArray(source, d_sourceArray, channelDesc));

	printf("Initiated Volume...\n");
}

extern "C" 
int runSimpleTexture(dim3 gridSize, dim3 blockSize, int * d_output,float * d_matrix,int imageX,int imageY,int imageZ,int modX, int modY)
{
	kernel_register<<<gridSize,blockSize>>>(d_output,d_matrix,imageX,imageY,imageZ,modX,modY);
	return 1;
}

extern "C" 
int runSlowTexture(dim3 gridSize, dim3 blockSize, int * d_output,float * d_matrix,int imageX,int imageY,int imageZ,int modX, int modY)
{
	kernel_atomicRegister<<<gridSize,blockSize>>>(d_output,d_matrix,imageX,imageY,imageZ,modX,modY);
	return 1;
}

extern "C" 
int runRegistrationWithSamples(dim3 gridSize, dim3 blockSize, int * d_output,int * d_samplesArray, float * d_matrix,int imageX,int imageY,int imageZ,int modX, int modY)
{
	kernel_registerSamples<<<gridSize,blockSize>>>(d_output,d_samplesArray,d_matrix,imageX,imageY,imageZ,modX,modY);
	return 1;
}

extern "C" 
int runOutputNormSourceVals(dim3 gridSize, dim3 blockSize, float * d_output,float * d_matrix,int imageX,int imageY,int imageZ,int modX, int modY)
{
	kernel_outputNormSourceVals<<<gridSize,blockSize>>>(d_output,d_matrix,imageX,imageY,imageZ,modX,modY);
	return 1;
}