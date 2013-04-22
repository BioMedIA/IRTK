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



#include "reduction/reduction.h"

template <class VoxelType> Reduction<VoxelType>::Reduction(VoxelType* hostArray, int size)
{
	temporaryBuffer = new DeviceArray<VoxelType>(hostArray, size);
	numInputs = size;
}

template <class VoxelType> Reduction<VoxelType>::~Reduction()
{

}

template <class VoxelType> 
void Reduction<VoxelType>::Sum(VoxelType * inputs, VoxelType * output, int numInputs, VoxelType multiplyFactor, cudaStream_t stream) {
	int blockSize = NumberThreadsPerBlockThatBestFit(numInputs, OPTIMAL_BLOCK_SIZE_REDUCTION);

	if (numInputs > SIZE_SMALL_CUDA_VECTOR) {
		int blocks = NumberBlocks(numInputs, blockSize);
		if (temporaryBuffer->Length() < blocks) temporaryBuffer->ResizeWithoutPreservingData(blocks);

		KernelSum<VoxelType>(stream, blocks, blockSize, inputs, temporaryBuffer->Pointer(), numInputs);

		inputs = temporaryBuffer->Pointer();
		numInputs = blocks;

		blockSize = NumberThreadsPerBlockThatBestFit(numInputs);
	}
			
	KernelSumSmallArray<VoxelType>(stream, blockSize, inputs, output, numInputs, multiplyFactor);
}

template <class VoxelType> 
void Reduction<VoxelType>::Min(VoxelType * inputs, VoxelType * output, int numInputs, cudaStream_t stream) {
	int blockSize = NumberThreadsPerBlockThatBestFit(numInputs, OPTIMAL_BLOCK_SIZE_REDUCTION);

	if (numInputs > SIZE_SMALL_CUDA_VECTOR) {
		int blocks = NumberBlocks(numInputs, blockSize);
		if (temporaryBuffer->Length() < blocks) temporaryBuffer->ResizeWithoutPreservingData(blocks);

		KernelMin<VoxelType>(stream, blocks, blockSize, inputs, temporaryBuffer->Pointer(), numInputs);

		inputs = temporaryBuffer->Pointer();
		numInputs = blocks;

		blockSize = NumberThreadsPerBlockThatBestFit(numInputs);
	}

	KernelMin<VoxelType>(stream, 1, blockSize, inputs, output, numInputs);
}

template <class VoxelType> 
void Reduction<VoxelType>::MinIndex(VoxelType * inputs, VoxelType * output, int * minIndex, int numInputs, cudaStream_t stream) {
	int blockSize = NumberThreadsPerBlockThatBestFit(numInputs, OPTIMAL_BLOCK_SIZE_REDUCTION);

	int * indexes = nullptr;

	if (numInputs > SIZE_SMALL_CUDA_VECTOR) {
		int blocks = NumberBlocks(numInputs, blockSize);

		int minSizeBuffer = blocks + (int) ceil(blocks * (sizeof(int) / (float) sizeof(float)));
		if (temporaryBuffer->Length() < minSizeBuffer) temporaryBuffer->ResizeWithoutPreservingData(minSizeBuffer);

		indexes = (int *)(temporaryBuffer->Pointer() + blocks);

		KernelMinIndexes<VoxelType>(stream, blocks, blockSize, inputs, temporaryBuffer->Pointer(), indexes, numInputs, nullptr);

		inputs = temporaryBuffer->Pointer();
		numInputs = blocks;

		blockSize = NumberThreadsPerBlockThatBestFit(numInputs);
	}

	KernelMinIndexes<VoxelType>(stream, 1, blockSize, inputs, output, minIndex, numInputs, indexes);
}

template <class VoxelType> 
void Reduction<VoxelType>::Max(VoxelType * inputs, VoxelType * output, int numInputs, cudaStream_t stream) {
	int blockSize = NumberThreadsPerBlockThatBestFit(numInputs, OPTIMAL_BLOCK_SIZE_REDUCTION);

	if (numInputs > SIZE_SMALL_CUDA_VECTOR) {
		int blocks = NumberBlocks(numInputs, blockSize);
		if (temporaryBuffer->Length() < blocks) temporaryBuffer->ResizeWithoutPreservingData(blocks);

		KernelMax<VoxelType>(stream, blocks, blockSize, inputs, temporaryBuffer->Pointer(), numInputs);

		inputs = temporaryBuffer->Pointer();
		numInputs = blocks;

		blockSize = NumberThreadsPerBlockThatBestFit(numInputs);
	}

	KernelMax<VoxelType>(stream, 1, blockSize, inputs, output, numInputs);
}

template <class VoxelType> 
void Reduction<VoxelType>::MaxIndex(VoxelType * inputs, VoxelType * output, int * maxIndex, int numInputs, cudaStream_t stream) {
	int blockSize = NumberThreadsPerBlockThatBestFit(numInputs, OPTIMAL_BLOCK_SIZE_REDUCTION);

	int * indexes = nullptr;

	if (numInputs > SIZE_SMALL_CUDA_VECTOR) {
		int blocks = NumberBlocks(numInputs, blockSize);

		int minSizeBuffer = blocks + (int) ceil(blocks * (sizeof(int) / (float) sizeof(float)));
		if (temporaryBuffer->Length() < minSizeBuffer) temporaryBuffer->ResizeWithoutPreservingData(minSizeBuffer);

		indexes = (int *)(temporaryBuffer->Pointer() + blocks);

		KernelMaxIndexes<VoxelType>(stream, blocks, blockSize, inputs, temporaryBuffer->Pointer(), indexes, numInputs, nullptr);

		inputs = temporaryBuffer->Pointer();
		numInputs = blocks;

		blockSize = NumberThreadsPerBlockThatBestFit(numInputs);
	}

	KernelMaxIndexes<VoxelType>(stream, 1, blockSize, inputs, output, maxIndex, numInputs, indexes);
}


template class irtkCULib_DLLAPI Reduction<char>;
template class irtkCULib_DLLAPI Reduction<unsigned char>;
template class irtkCULib_DLLAPI Reduction<short>;
template class irtkCULib_DLLAPI Reduction<unsigned short>;
template class irtkCULib_DLLAPI Reduction<int>;
template class irtkCULib_DLLAPI Reduction<unsigned int>;
template class irtkCULib_DLLAPI Reduction<float>;
template class irtkCULib_DLLAPI Reduction<double>;