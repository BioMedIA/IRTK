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


#include "reduction/reduction.h"
#include <limits>

template <class VoxelType, int blockSize> KERNEL Max(VoxelType * inputs, VoxelType * output, int numInputs, VoxelType typeMin) {
	//extern __shared__ VoxelType maxvalue[];
	VoxelType *maxvalue = SharedMemory<VoxelType>();

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	VoxelType value = typeMin;
	if (idx < numInputs) value = inputs[idx];

	maxvalue[threadIdx.x] = value;
	__syncthreads();

	if (blockSize >= 1024) {
		if (threadIdx.x < 512 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 512]) maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 512];
		__syncthreads();
	}

	if (blockSize >= 512) {
		if (threadIdx.x < 256 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 256]) maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 256];
		__syncthreads();
	}

	if (blockSize >= 256) {
		if (threadIdx.x < 128 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 128]) maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 128];
		__syncthreads();
	}

	if (blockSize >= 128) {
		if (threadIdx.x < 64 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 64]) maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 64];
		__syncthreads();
	}

	if (threadIdx.x < 32) {
		volatile VoxelType * _maxvalue = maxvalue;

		if (blockSize >= 64) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 32]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 32];
		}

		if (blockSize >= 32) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 16]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 16];
		}

		if (blockSize >= 16) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 8]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 8];
		}

		if (blockSize >= 8) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 4]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 4];
		}

		if (blockSize >= 4) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 2]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 2];
		}

		if (blockSize >= 2) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 1]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 1];
		}

		if (threadIdx.x == 0) {
			output[blockIdx.x] = maxvalue[0];
		}
	}
}

template <class VoxelType, int blockSize> KERNEL MaxIndex(VoxelType * inputs, VoxelType * output, int * indexes, int numInputs, VoxelType typeMin) {
	//extern __shared__ VoxelType maxvalue[];
	VoxelType *maxvalue = SharedMemory<VoxelType>();

	int * maxpos = (int *) (maxvalue + blockDim.x);

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	VoxelType value = typeMin;
	if (idx < numInputs) value = inputs[idx];

	maxvalue[threadIdx.x] = value;
	maxpos[threadIdx.x] = idx;
	__syncthreads();

	if (blockSize >= 1024) {
		if (threadIdx.x < 512 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 512]) {
			maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 512];
			maxpos[threadIdx.x] = maxpos[threadIdx.x + 512];
		}
		__syncthreads();
	}

	if (blockSize >= 512) {
		if (threadIdx.x < 256 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 256]) {
			maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 256];
			maxpos[threadIdx.x] = maxpos[threadIdx.x + 256];
		}
		__syncthreads();
	}

	if (blockSize >= 256) {
		if (threadIdx.x < 128 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 128]) {
			maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 128];
			maxpos[threadIdx.x] = maxpos[threadIdx.x + 128];
		}
		__syncthreads();
	}

	if (blockSize >= 128) {
		if (threadIdx.x < 64 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 64]) {
			maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 64];
			maxpos[threadIdx.x] = maxpos[threadIdx.x + 64];
		}
		__syncthreads();
	}

	if (threadIdx.x < 32) {
		volatile VoxelType * _maxvalue = maxvalue;
		volatile int * _maxpos = maxpos;

		if (blockSize >= 64) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 32]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 32];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 32];
			}
		}

		if (blockSize >= 32) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 16]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 16];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 16];
			}
		}

		if (blockSize >= 16) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 8]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 8];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 8];
			}
		}

		if (blockSize >= 8) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 4]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 4];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 4];
			}
		}

		if (blockSize >= 4) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 2]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 2];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 2];
			}
		}

		if (blockSize >= 2) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 1]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 1];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 1];
			}
		}

		if (threadIdx.x == 0) {
			output[blockIdx.x] = maxvalue[0];
			indexes[blockIdx.x] = maxpos[0];
		}
	}
}

template <class VoxelType, int blockSize> KERNEL MaxSmallArray(VoxelType * inputs, VoxelType * output, int numInputs, VoxelType typeMin) {
	//extern __shared__ VoxelType maxvalue[];
	VoxelType *maxvalue = SharedMemory<VoxelType>();

	maxvalue[threadIdx.x] = typeMin;
	for(int i = threadIdx.x; i < numInputs; i += blockDim.x) if (maxvalue[threadIdx.x] < inputs[i]) maxvalue[threadIdx.x] = inputs[i];
	__syncthreads();

	if (blockSize >= 1024) {
		if (threadIdx.x < 512 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 512]) maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 512];
		__syncthreads();
	}

	if (blockSize >= 512) {
		if (threadIdx.x < 256 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 256]) maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 256];
		__syncthreads();
	}

	if (blockSize >= 256) {
		if (threadIdx.x < 128 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 128]) maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 128];
		__syncthreads();
	}

	if (blockSize >= 128) {
		if (threadIdx.x < 64 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 64]) maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 64];
		__syncthreads();
	}

	if (threadIdx.x < 32) {
		volatile VoxelType * _maxvalue = maxvalue;

		if (blockSize >= 64) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 32]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 32];
		}

		if (blockSize >= 32) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 16]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 16];
		}

		if (blockSize >= 16) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 8]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 8];
		}

		if (blockSize >= 8) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 4]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 4];
		}

		if (blockSize >= 4) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 2]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 2];
		}

		if (blockSize >= 2) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 1]) _maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 1];
		}

		if (threadIdx.x == 0) {
			output[blockIdx.x] = maxvalue[0];
		}
	}
}

template <class VoxelType, int blockSize> KERNEL MaxSmallArrayIndex(VoxelType * inputs, VoxelType * output, int * maxIndex, int numInputs, int * indexes, VoxelType typeMin) {
	//extern __shared__ VoxelType maxvalue[];
	VoxelType *maxvalue = SharedMemory<VoxelType>();

	int * maxpos = (int *) (maxvalue + blockDim.x);

	maxvalue[threadIdx.x] = typeMin;
	for(int i = threadIdx.x; i < numInputs; i += blockDim.x) {
		if (maxvalue[threadIdx.x] < inputs[i]) {
			maxvalue[threadIdx.x] = inputs[i];
			if (indexes != nullptr) {
				maxpos[threadIdx.x] = indexes[i];
			} else {
				maxpos[threadIdx.x] = i;
			}
		}
	}
	__syncthreads();

	if (blockSize >= 1024) {
		if (threadIdx.x < 512 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 512]) {
			maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 512];
			maxpos[threadIdx.x] = maxpos[threadIdx.x + 512];
		}
		__syncthreads();
	}

	if (blockSize >= 512) {
		if (threadIdx.x < 256 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 256]) {
			maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 256];
			maxpos[threadIdx.x] = maxpos[threadIdx.x + 256];
		}
		__syncthreads();
	}

	if (blockSize >= 256) {
		if (threadIdx.x < 128 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 128]) {
			maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 128];
			maxpos[threadIdx.x] = maxpos[threadIdx.x + 128];
		}
		__syncthreads();
	}

	if (blockSize >= 128) {
		if (threadIdx.x < 64 && maxvalue[threadIdx.x] < maxvalue[threadIdx.x + 64]) {
			maxvalue[threadIdx.x] = maxvalue[threadIdx.x + 64];
			maxpos[threadIdx.x] = maxpos[threadIdx.x + 64];
		}
		__syncthreads();
	}

	if (threadIdx.x < 32) {
		volatile VoxelType * _maxvalue = maxvalue;
		volatile int * _maxpos = maxpos;

		if (blockSize >= 64) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 32]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 32];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 32];
			}
		}

		if (blockSize >= 32) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 16]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 16];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 16];
			}
		}

		if (blockSize >= 16) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 8]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 8];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 8];
			}
		}

		if (blockSize >= 8) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 4]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 4];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 4];
			}
		}

		if (blockSize >= 4) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 2]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 2];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 2];
			}
		}

		if (blockSize >= 2) {
			if (_maxvalue[threadIdx.x] < _maxvalue[threadIdx.x + 1]) {
				_maxvalue[threadIdx.x] = _maxvalue[threadIdx.x + 1];
				_maxpos[threadIdx.x] = _maxpos[threadIdx.x + 1];
			}
		}

		if (threadIdx.x == 0) {
			output[blockIdx.x] = maxvalue[0];
			maxIndex[blockIdx.x] = maxpos[0];
		}
	}
}

template <class VoxelType> 
void KernelMax(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * output, int numInputs) {
	VoxelType typeMin = numeric_limits<VoxelType>::min();
	if (blocks == 1) {
		switch(blockSize) {
#ifdef FERMI
		case 1024:
			MaxSmallArray<VoxelType, 1024><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
#endif
		case 512:
			MaxSmallArray<VoxelType, 512><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 256:
			MaxSmallArray<VoxelType, 256><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 128:
			MaxSmallArray<VoxelType, 128><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 64:
			MaxSmallArray<VoxelType, 64><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 32:
			MaxSmallArray<VoxelType, 32><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 16:
			MaxSmallArray<VoxelType, 16><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 8:
			MaxSmallArray<VoxelType, 8><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 4:
			MaxSmallArray<VoxelType, 4><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 2:
			MaxSmallArray<VoxelType, 2><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 1:
			MaxSmallArray<VoxelType, 1><<<1, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		}
	} else {
		switch(blockSize) {
#ifdef FERMI
		case 1024:
			Max<VoxelType, 1024><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
#endif
		case 512:
			Max<VoxelType, 512><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 256:
			Max<VoxelType, 256><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 128:
			Max<VoxelType, 128><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 64:
			Max<VoxelType, 64><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 32:
			Max<VoxelType, 32><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 16:
			Max<VoxelType, 16><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 8:
			Max<VoxelType, 8><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 4:
			Max<VoxelType, 4><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 2:
			Max<VoxelType, 2><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		case 1:
			Max<VoxelType, 1><<<blocks, blockSize, blockSize * sizeof(VoxelType), stream>>>(inputs, output, numInputs, typeMin);
			break;
		}
	}
}

template <class VoxelType> 
void KernelMaxIndexes(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * output, int * maxIndexes, int numInputs, int * indexes) {
	VoxelType typeMin = numeric_limits<VoxelType>::min();
	if (blocks == 1) {
		switch(blockSize) {
#ifdef FERMI
		case 1024:
			MaxSmallArrayIndex<VoxelType, 1024><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
#endif
		case 512:
			MaxSmallArrayIndex<VoxelType, 512><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		case 256:
			MaxSmallArrayIndex<VoxelType, 256><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		case 128:
			MaxSmallArrayIndex<VoxelType, 128><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		case 64:
			MaxSmallArrayIndex<VoxelType, 64><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		case 32:
			MaxSmallArrayIndex<VoxelType, 32><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		case 16:
			MaxSmallArrayIndex<VoxelType, 16><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		case 8:
			MaxSmallArrayIndex<VoxelType, 8><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		case 4:
			MaxSmallArrayIndex<VoxelType, 4><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		case 2:
			MaxSmallArrayIndex<VoxelType, 2><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		case 1:
			MaxSmallArrayIndex<VoxelType, 1><<<1, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, indexes, typeMin);
			break;
		}
	} else {
		switch(blockSize) {
#ifdef FERMI
		case 1024:
			MaxIndex<VoxelType, 1024><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
#endif
		case 512:
			MaxIndex<VoxelType, 512><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		case 256:
			MaxIndex<VoxelType, 256><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		case 128:
			MaxIndex<VoxelType, 128><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		case 64:
			MaxIndex<VoxelType, 64><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		case 32:
			MaxIndex<VoxelType, 32><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		case 16:
			MaxIndex<VoxelType, 16><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		case 8:
			MaxIndex<VoxelType, 8><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		case 4:
			MaxIndex<VoxelType, 4><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		case 2:
			MaxIndex<VoxelType, 2><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		case 1:
			MaxIndex<VoxelType, 1><<<blocks, blockSize, blockSize * (sizeof(VoxelType) + sizeof(int)), stream>>>(inputs, output, maxIndexes, numInputs, typeMin);
			break;
		}
	}
}


#define INSTANCE_MACRO(TYPE) template void KernelMaxIndexes<TYPE>(cudaStream_t stream, int blocks, int blockSize, TYPE * inputs, TYPE * output, int * minIndexes, int numInputs, int * indexes);
INSTANCE_MACRO(char)
	INSTANCE_MACRO(unsigned char)
	INSTANCE_MACRO(short)
	INSTANCE_MACRO(unsigned short)
	INSTANCE_MACRO(int)
	INSTANCE_MACRO(unsigned int)
	INSTANCE_MACRO(float)
	INSTANCE_MACRO(double)

#define INSTANCE_MACRO1(TYPE) template void KernelMax<TYPE>(cudaStream_t stream, int blocks, int blockSize, TYPE * inputs, TYPE * output, int numInputs);
	INSTANCE_MACRO1(char)
	INSTANCE_MACRO1(unsigned char)
	INSTANCE_MACRO1(short)
	INSTANCE_MACRO1(unsigned short)
	INSTANCE_MACRO1(int)
	INSTANCE_MACRO1(unsigned int)
	INSTANCE_MACRO1(float)
	INSTANCE_MACRO1(double)