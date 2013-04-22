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



#include "reduction/SumWarp.h"
#include "reduction/reduction.h"

template <class VoxelType, int blockSize> KERNEL Sum(VoxelType * inputs, VoxelType * outputs, int numInputs) {

	//extern __shared__ VoxelType sum[];
	VoxelType *sum = SharedMemory<VoxelType>();

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	VoxelType value = 0.0;
	if (idx < numInputs) value = inputs[idx];

	sum[threadIdx.x] = value;
	__syncthreads();

	SumBeforeWarp<VoxelType, blockSize>(sum);

	if (threadIdx.x < 32) {
		SumWarp<VoxelType, blockSize>(sum);
		if (threadIdx.x == 0) outputs[blockIdx.x] = sum[0];
	}
}

template <class VoxelType> 
void KernelSum(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * outputs, int numInputs) {
	switch(blockSize) {
#ifdef FERMI
	case 1024:
		Sum<VoxelType, 1024><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
#endif
	case 512:
		Sum<VoxelType, 512><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	case 256:
		Sum<VoxelType, 256><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	case 128:
		Sum<VoxelType, 128><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	case 64:
		Sum<VoxelType, 64><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	case 32:
		Sum<VoxelType, 32><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	case 16:
		Sum<VoxelType, 16><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	case 8:
		Sum<VoxelType, 8><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	case 4:
		Sum<VoxelType, 4><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	case 2:
		Sum<VoxelType, 2><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	case 1:
		Sum<VoxelType, 1><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, outputs, numInputs);
		break;
	}
}

template <class VoxelType, int blockSize> KERNEL SumSmallArray(VoxelType * inputs, VoxelType * output, int numInputs, VoxelType multiplyFactor) {
	//extern __shared__ VoxelType sum[];
	VoxelType *sum = SharedMemory<VoxelType>();

	VoxelType value = 0.0;
	for(int i = threadIdx.x; i < numInputs; i += blockDim.x) value += inputs[i]; 
	sum[threadIdx.x] = value;
	__syncthreads();

	SumBeforeWarp<VoxelType, blockSize>(sum);

	if (threadIdx.x < 32) {
		SumWarp<VoxelType, blockSize>(sum);

		if (threadIdx.x == 0) output[blockIdx.x] = sum[0] * multiplyFactor;
	}
}

template <class VoxelType> 
void KernelSumSmallArray(cudaStream_t stream, int blockSize, VoxelType * inputs, VoxelType * output, int numInputs, VoxelType multiplyFactor) {
	switch(blockSize) {
#ifdef FERMI
	case 1024:
		SumSmallArray<VoxelType, 1024><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
#endif
	case 512:
		SumSmallArray<VoxelType,512><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	case 256:
		SumSmallArray<VoxelType,256><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	case 128:
		SumSmallArray<VoxelType,128><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	case 64:
		SumSmallArray<VoxelType,64><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	case 32:
		SumSmallArray<VoxelType,32><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	case 16:
		SumSmallArray<VoxelType,16><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	case 8:
		SumSmallArray<VoxelType,8><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	case 4:
		SumSmallArray<VoxelType,4><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	case 2:
		SumSmallArray<VoxelType,2><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	case 1:
		SumSmallArray<VoxelType,1><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, multiplyFactor);
		break;
	}
}

#define INSTANCE_MACRO(TYPE) template void KernelSumSmallArray<TYPE>(cudaStream_t stream, int blockSize, TYPE * inputs, TYPE * output, int numInputs, TYPE multiplyFactor);
INSTANCE_MACRO(char)
	INSTANCE_MACRO(unsigned char)
	INSTANCE_MACRO(short)
	INSTANCE_MACRO(unsigned short)
	INSTANCE_MACRO(int)
	INSTANCE_MACRO(unsigned int)
	INSTANCE_MACRO(float)
	INSTANCE_MACRO(double)

#define INSTANCE_MACRO1(TYPE) template void KernelSum<TYPE>(cudaStream_t stream, int blocks, int blockSize, TYPE * inputs, TYPE * outputs, int numInputs);
	INSTANCE_MACRO1(char)
	INSTANCE_MACRO1(unsigned char)
	INSTANCE_MACRO1(short)
	INSTANCE_MACRO1(unsigned short)
	INSTANCE_MACRO1(int)
	INSTANCE_MACRO1(unsigned int)
	INSTANCE_MACRO1(float)
	INSTANCE_MACRO1(double)