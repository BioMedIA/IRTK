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


#include <limits>
#include "reduction/reduction.h"

template <class VoxelType, int blockSize> KERNEL Min(VoxelType * inputs, VoxelType * output, int numInputs, VoxelType typeMax) {
	//extern __shared__ VoxelType minvalue[];
	VoxelType *minvalue = SharedMemory<VoxelType>();

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	VoxelType value = typeMax;
	if (idx < numInputs) value = inputs[idx];

	minvalue[threadIdx.x] = value;
	__syncthreads();

	if (blockSize >= 1024) {
		if (threadIdx.x < 512 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 512]) minvalue[threadIdx.x] = minvalue[threadIdx.x + 512];
		__syncthreads();
	}

	if (blockSize >= 512) {
		if (threadIdx.x < 256 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 256]) minvalue[threadIdx.x] = minvalue[threadIdx.x + 256];
		__syncthreads();
	}

	if (blockSize >= 256) {
		if (threadIdx.x < 128 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 128]) minvalue[threadIdx.x] = minvalue[threadIdx.x + 128];
		__syncthreads();
	}

	if (blockSize >= 128) {
		if (threadIdx.x < 64 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 64]) minvalue[threadIdx.x] = minvalue[threadIdx.x + 64];
		__syncthreads();
	}

	if (threadIdx.x < 32) {
		volatile VoxelType * _minvalue = minvalue;

		if (blockSize >= 64) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 32]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 32];
		}

		if (blockSize >= 32) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 16]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 16];
		}

		if (blockSize >= 16) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 8]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 8];
		}

		if (blockSize >= 8) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 4]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 4];
		}

		if (blockSize >= 4) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 2]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 2];
		}

		if (blockSize >= 2) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 1]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 1];
		}

		if (threadIdx.x == 0) {
			output[blockIdx.x] = minvalue[0];
		}
	}
}

template <class VoxelType, int blockSize> KERNEL MinIndex(VoxelType * inputs, VoxelType * output, int * indexes, int numInputs, VoxelType typeMax) {
	//extern __shared__ VoxelType minvalue[];
	VoxelType *minvalue = SharedMemory<VoxelType>();

	int * minpos = (int *) (minvalue + blockDim.x);

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	VoxelType value = typeMax;
	if (idx < numInputs) value = inputs[idx];

	minvalue[threadIdx.x] = value;
	minpos[threadIdx.x] = idx;
	__syncthreads();

	if (blockSize >= 1024) {
		if (threadIdx.x < 512 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 512]) {
			minvalue[threadIdx.x] = minvalue[threadIdx.x + 512];
			minpos[threadIdx.x] = minpos[threadIdx.x + 512];
		}
		__syncthreads();
	}

	if (blockSize >= 512) {
		if (threadIdx.x < 256 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 256]) {
			minvalue[threadIdx.x] = minvalue[threadIdx.x + 256];
			minpos[threadIdx.x] = minpos[threadIdx.x + 256];
		}
		__syncthreads();
	}

	if (blockSize >= 256) {
		if (threadIdx.x < 128 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 128]) {
			minvalue[threadIdx.x] = minvalue[threadIdx.x + 128];
			minpos[threadIdx.x] = minpos[threadIdx.x + 128];
		}
		__syncthreads();
	}

	if (blockSize >= 128) {
		if (threadIdx.x < 64 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 64]) {
			minvalue[threadIdx.x] = minvalue[threadIdx.x + 64];
			minpos[threadIdx.x] = minpos[threadIdx.x + 64];
		}
		__syncthreads();
	}

	if (threadIdx.x < 32) {
		volatile VoxelType * _minvalue = minvalue;
		volatile int * _minpos = minpos;

		if (blockSize >= 64) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 32]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 32];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 32];
			}
		}

		if (blockSize >= 32) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 16]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 16];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 16];
			}
		}

		if (blockSize >= 16) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 8]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 8];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 8];
			}
		}

		if (blockSize >= 8) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 4]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 4];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 4];
			}
		}

		if (blockSize >= 4) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 2]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 2];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 2];
			}
		}

		if (blockSize >= 2) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 1]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 1];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 1];
			}
		}

		if (threadIdx.x == 0) {
			output[blockIdx.x] = minvalue[0];
			indexes[blockIdx.x] = minpos[0];
		}
	}
}

template <class VoxelType, int blockSize> KERNEL MinSmallArray(VoxelType * inputs, VoxelType * output, int numInputs, VoxelType typeMax) {
	//extern __shared__ VoxelType minvalue[];
	VoxelType *minvalue = SharedMemory<VoxelType>();

	minvalue[threadIdx.x] = typeMax;
	for(int i = threadIdx.x; i < numInputs; i += blockDim.x) if (minvalue[threadIdx.x] > inputs[i]) minvalue[threadIdx.x] = inputs[i];
	__syncthreads();

	if (blockSize >= 1024) {
		if (threadIdx.x < 512 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 512]) minvalue[threadIdx.x] = minvalue[threadIdx.x + 512];
		__syncthreads();
	}

	if (blockSize >= 512) {
		if (threadIdx.x < 256 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 256]) minvalue[threadIdx.x] = minvalue[threadIdx.x + 256];
		__syncthreads();
	}

	if (blockSize >= 256) {
		if (threadIdx.x < 128 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 128]) minvalue[threadIdx.x] = minvalue[threadIdx.x + 128];
		__syncthreads();
	}

	if (blockSize >= 128) {
		if (threadIdx.x < 64 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 64]) minvalue[threadIdx.x] = minvalue[threadIdx.x + 64];
		__syncthreads();
	}

	if (threadIdx.x < 32) {
		volatile VoxelType * _minvalue = minvalue;

		if (blockSize >= 64) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 32]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 32];
		}

		if (blockSize >= 32) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 16]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 16];
		}

		if (blockSize >= 16) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 8]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 8];
		}

		if (blockSize >= 8) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 4]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 4];
		}

		if (blockSize >= 4) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 2]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 2];
		}

		if (blockSize >= 2) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 1]) _minvalue[threadIdx.x] = _minvalue[threadIdx.x + 1];
		}

		if (threadIdx.x == 0) {
			output[blockIdx.x] = minvalue[0];
		}
	}
}

template <class VoxelType, int blockSize> KERNEL MinSmallArrayIndex(VoxelType * inputs, VoxelType * output, int * minIndex, int numInputs, int * indexes, VoxelType typeMax) {
	//extern __shared__ VoxelType minvalue[];
	VoxelType *minvalue = SharedMemory<VoxelType>();

	int * minpos = (int *) (minvalue + blockDim.x);

	minvalue[threadIdx.x] = typeMax;
	for(int i = threadIdx.x; i < numInputs; i += blockDim.x) {
		if (minvalue[threadIdx.x] > inputs[i]) {
			minvalue[threadIdx.x] = inputs[i];
			if (indexes != nullptr) {
				minpos[threadIdx.x] = indexes[i];
			} else {
				minpos[threadIdx.x] = i;
			}
		}
	}
	__syncthreads();

	if (blockSize >= 1024) {
		if (threadIdx.x < 512 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 512]) {
			minvalue[threadIdx.x] = minvalue[threadIdx.x + 512];
			minpos[threadIdx.x] = minpos[threadIdx.x + 512];
		}
		__syncthreads();
	}

	if (blockSize >= 512) {
		if (threadIdx.x < 256 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 256]) {
			minvalue[threadIdx.x] = minvalue[threadIdx.x + 256];
			minpos[threadIdx.x] = minpos[threadIdx.x + 256];
		}
		__syncthreads();
	}

	if (blockSize >= 256) {
		if (threadIdx.x < 128 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 128]) {
			minvalue[threadIdx.x] = minvalue[threadIdx.x + 128];
			minpos[threadIdx.x] = minpos[threadIdx.x + 128];
		}
		__syncthreads();
	}

	if (blockSize >= 128) {
		if (threadIdx.x < 64 && minvalue[threadIdx.x] > minvalue[threadIdx.x + 64]) {
			minvalue[threadIdx.x] = minvalue[threadIdx.x + 64];
			minpos[threadIdx.x] = minpos[threadIdx.x + 64];
		}
		__syncthreads();
	}

	if (threadIdx.x < 32) {
		volatile VoxelType * _minvalue = minvalue;
		volatile int * _minpos = minpos;

		if (blockSize >= 64) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 32]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 32];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 32];
			}
		}

		if (blockSize >= 32) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 16]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 16];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 16];
			}
		}

		if (blockSize >= 16) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 8]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 8];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 8];
			}
		}

		if (blockSize >= 8) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 4]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 4];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 4];
			}
		}

		if (blockSize >= 4) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 2]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 2];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 2];
			}
		}

		if (blockSize >= 2) {
			if (_minvalue[threadIdx.x] > _minvalue[threadIdx.x + 1]) {
				_minvalue[threadIdx.x] = _minvalue[threadIdx.x + 1];
				_minpos[threadIdx.x] = _minpos[threadIdx.x + 1];
			}
		}

		if (threadIdx.x == 0) {
			output[blockIdx.x] = minvalue[0];
			minIndex[blockIdx.x] = minpos[0];
		}
	}
}

template <class VoxelType> 
void KernelMin(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * output, int numInputs) {
	VoxelType typeMax = numeric_limits<VoxelType>::max();
	if (blocks == 1) {
		switch(blockSize) {
#ifdef FERMI
		case 1024:
			MinSmallArray<VoxelType, 1024><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
#endif
		case 512:
			MinSmallArray<VoxelType, 512><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 256:
			MinSmallArray<VoxelType, 256><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 128:
			MinSmallArray<VoxelType, 128><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 64:
			MinSmallArray<VoxelType, 64><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 32:
			MinSmallArray<VoxelType, 32><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 16:
			MinSmallArray<VoxelType, 16><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 8:
			MinSmallArray<VoxelType, 8><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 4:
			MinSmallArray<VoxelType, 4><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 2:
			MinSmallArray<VoxelType, 2><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 1:
			MinSmallArray<VoxelType, 1><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		}
	} else {
		switch(blockSize) {
#ifdef FERMI
		case 1024:
			Min<VoxelType, 1024><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
#endif
		case 512:
			Min<VoxelType, 512><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 256:
			Min<VoxelType, 256><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 128:
			Min<VoxelType, 128><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 64:
			Min<VoxelType, 64><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 32:
			Min<VoxelType, 32><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 16:
			Min<VoxelType, 16><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 8:
			Min<VoxelType, 8><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 4:
			Min<VoxelType, 4><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 2:
			Min<VoxelType, 2><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		case 1:
			Min<VoxelType, 1><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMax);
			break;
		}
	}
}

template <class VoxelType> 
void KernelMinIndexes(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * output, int * minIndexes, int numInputs, int * indexes) {
	VoxelType typeMax = numeric_limits<VoxelType>::max();
	if (blocks == 1) {
		switch(blockSize) {
#ifdef FERMI
		case 1024:
			MinSmallArrayIndex<VoxelType, 1024><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
#endif
		case 512:
			MinSmallArrayIndex<VoxelType, 512><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		case 256:
			MinSmallArrayIndex<VoxelType, 256><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		case 128:
			MinSmallArrayIndex<VoxelType, 128><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		case 64:
			MinSmallArrayIndex<VoxelType, 64><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		case 32:
			MinSmallArrayIndex<VoxelType, 32><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		case 16:
			MinSmallArrayIndex<VoxelType, 16><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		case 8:
			MinSmallArrayIndex<VoxelType, 8><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		case 4:
			MinSmallArrayIndex<VoxelType, 4><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		case 2:
			MinSmallArrayIndex<VoxelType, 2><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		case 1:
			MinSmallArrayIndex<VoxelType, 1><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, indexes, typeMax);
			break;
		}
	} else {
		switch(blockSize) {
#ifdef FERMI
		case 1024:
			MinIndex<VoxelType, 1024><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
#endif
		case 512:
			MinIndex<VoxelType, 512><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		case 256:
			MinIndex<VoxelType, 256><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		case 128:
			MinIndex<VoxelType, 128><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		case 64:
			MinIndex<VoxelType, 64><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		case 32:
			MinIndex<VoxelType, 32><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		case 16:
			MinIndex<VoxelType, 16><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		case 8:
			MinIndex<VoxelType, 8><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		case 4:
			MinIndex<VoxelType, 4><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		case 2:
			MinIndex<VoxelType, 2><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		case 1:
			MinIndex<VoxelType, 1><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, minIndexes, numInputs, typeMax);
			break;
		}
	}
}


#define INSTANCE_MACRO(TYPE) template void KernelMinIndexes<TYPE>(cudaStream_t stream, int blocks, int blockSize, TYPE * inputs, TYPE * output, int * minIndexes, int numInputs, int * indexes);
INSTANCE_MACRO(char)
	INSTANCE_MACRO(unsigned char)
	INSTANCE_MACRO(short)
	INSTANCE_MACRO(unsigned short)
	INSTANCE_MACRO(int)
	INSTANCE_MACRO(unsigned int)
	INSTANCE_MACRO(float)
	INSTANCE_MACRO(double)

#define INSTANCE_MACRO1(TYPE) template void KernelMin<TYPE>(cudaStream_t stream, int blocks, int blockSize, TYPE * inputs, TYPE * output, int numInputs);
	INSTANCE_MACRO1(char)
	INSTANCE_MACRO1(unsigned char)
	INSTANCE_MACRO1(short)
	INSTANCE_MACRO1(unsigned short)
	INSTANCE_MACRO1(int)
	INSTANCE_MACRO1(unsigned int)
	INSTANCE_MACRO1(float)
	INSTANCE_MACRO1(double)