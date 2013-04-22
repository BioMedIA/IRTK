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


#ifndef GPUMLib_reduction_h
#define GPUMLib_reduction_h

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>
#include <cmath>

#include "CudaDefinitions.h"
#include "Utilities.h"
#include "../memory/CudaArray.h"
#include "../memory/DeviceMatrix.h"
#include "../memory/DeviceAccessibleVariable.h"

#include "irtkCUSharedLibMacros.h"

using namespace std;


// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
#ifdef __CUDACC__
template<class T>
struct SharedMemory
{
	__device__ inline operator       T *()
	{
		extern __shared__ int __smem[];
		return (T *)__smem;
	}

	__device__ inline operator const T *() const
	{
		extern __shared__ int __smem[];
		return (T *)__smem;
	}
};

// specialize for double to avoid unaligned memory
// access compile errors
template<>
struct SharedMemory<double>
{
	__device__ inline operator       double *()
	{
		extern __shared__ double __smem_d[];
		return (double *)__smem_d;
	}

	__device__ inline operator const double *() const
	{
		extern __shared__ double __smem_d[];
		return (double *)__smem_d;
	}
};
#endif

//! \addtogroup reduction Reduction framework
//! @{

//! Kernel to sum an array. For small arrays use KernelSumSmallArray instead.
//! \param[in] stream CUDA stream
//! \param[in] blocks Number of thread blocks 
//! \param[in] blockSize Block size (number of threads per block)
//! \param[in] inputs Values to be summed
//! \param[out] outputs Array that will contain the partial sums of each block
//! \param[in] numInputs Number of inputs
//! \sa KernelSumSmallArray, SIZE_SMALL_CUDA_VECTOR
template <class VoxelType>
void KernelSum(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * outputs, int numInputs);

//! Kernel to sum a small array, multiply the result by a given factor and place the result in the output.
//! \param[in] stream CUDA stream
//! \param[in] blockSize Block size (number of threads per block)
//! \param[in] inputs Values to be summed
//! \param[out] output Pointer to the location that will contain the sum output
//! \param[in] numInputs Number of inputs
//! \param[in] multiplyFactor Multiply factor (optional, by default 1.0)
//! \sa KernelSum, SIZE_SMALL_CUDA_VECTOR
template <class VoxelType>
void KernelSumSmallArray(cudaStream_t stream, int blockSize, VoxelType * inputs, VoxelType * output, int numInputs, VoxelType multiplyFactor);

//! Kernel to compute the minimum of an array. 
//! \param[in] stream CUDA stream
//! \param[in] blocks Number of thread blocks
//! \param[in] blockSize Block size (number of threads per block)
//! \param[in] inputs input array
//! \param[out] output Pointer to the location that will contain the minimum
//! \param[in] numInputs Number of inputs
template <class VoxelType>
void KernelMin(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * output, int numInputs);

//! Kernel to compute the minimum of an array and its index within the array. 
//! \param[in] stream CUDA stream
//! \param[in] blocks Number of thread blocks
//! \param[in] blockSize Block size (number of threads per block)
//! \param[in] inputs input array
//! \param[out] output Pointer to the location that will contain the minimum
//! \param[out] minIndexes Pointer to the location that will contain the index of one of the minimums
//! \param[in] numInputs Number of inputs
//! \param[in] indexes Buffer used to tempory store the indexes. Must have the same size of the inputs array.
template <class VoxelType>
void KernelMinIndexes(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * output, int * minIndexes, int numInputs, int * indexes);

//! Kernel to compute the maximum of an array. 
//! \param[in] stream CUDA stream
//! \param[in] blocks Number of thread blocks
//! \param[in] blockSize Block size (number of threads per block)
//! \param[in] inputs input array
//! \param[out] output Pointer to the location that will contain the maximum
//! \param[in] numInputs Number of inputs
template <class VoxelType>
void KernelMax(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * output, int numInputs);

//! Kernel to compute the maximum of an array and its index within the array. 
//! \param[in] stream CUDA stream
//! \param[in] blocks Number of thread blocks
//! \param[in] blockSize Block size (number of threads per block)
//! \param[in] inputs input array
//! \param[out] output Pointer to the location that will contain the maximum
//! \param[out] maxIndexes Pointer to the location that will contain the index of one of the maximums
//! \param[in] numInputs Number of inputs
//! \param[in] indexes Buffer used to tempory store the indexes. Must have the same size of the inputs array.
template <class VoxelType>
void KernelMaxIndexes(cudaStream_t stream, int blocks, int blockSize, VoxelType * inputs, VoxelType * output, int * maxIndexes, int numInputs, int * indexes);

//! Provides reduction functions (Sum, Average, Max, Min, ...).
template <class VoxelType>
class Reduction {
	private:
		void Sum(VoxelType * inputs, VoxelType * output, int numInputs, VoxelType multiplyFactor, cudaStream_t stream);

		void MinIndex(VoxelType * inputs, VoxelType * output, int * minIndex, int numInputs, cudaStream_t stream);
		void Min(VoxelType * inputs, VoxelType * output, int numInputs, cudaStream_t stream);

		void Max(VoxelType * inputs, VoxelType * output, int numInputs, cudaStream_t stream);
		void MaxIndex(VoxelType * inputs, VoxelType * output, int * minIndex, int numInputs, cudaStream_t stream);

		int numInputs;

	public:
		
		//Constructor
		Reduction(VoxelType* hostArray, int size);

		//Destructor
		~Reduction();

		//! Temporary buffer used for the reduction tasks. Programmers may take advantage of it for other tasks (hence, it is declared as public).
		DeviceArray<VoxelType>* temporaryBuffer;

		//! Sums all the elements of an input array, multiplies the sum by a given factor and places the result in the output
		//! \param[in] inputs Values to be summed
		//! \param[out] output Pointer to the memory address that will contain the sum output
		//! \param[in] multiplyFactor Multiply factor (optional, by default 1.0)
		//! \param[in] stream CUDA stream (optional)
		//void Sum(VoxelType * output, VoxelType multiplyFactor = 1.0, cudaStream_t stream = nullptr) {
		//	Sum(temporaryBuffer->Pointer(), output, temporaryBuffer->Length(), multiplyFactor, stream);
		//}

		//! Sums all the elements of an input array, multiplies the sum by a given factor and places the result in the output
		//! \param[in] inputs Values to be summed
		//! \param[out] output Array that will contain the sum output (in position 0)
		//! \param[in] multiplyFactor Multiply factor (optional, by default 1.0)
		//! \param[in] stream CUDA stream (optional)
		void Sum(VoxelType * output, VoxelType multiplyFactor = 1.0, cudaStream_t stream = nullptr) {
			DeviceArray<VoxelType> d_out(*temporaryBuffer);
			Sum(temporaryBuffer->Pointer(), d_out.Pointer(), temporaryBuffer->Length(), multiplyFactor, stream);
			checkCudaErrors(cudaMemcpy(output,d_out.Pointer(),sizeof(VoxelType),cudaMemcpyDeviceToHost));
		}

		//! Averages the elements of an input array, placing the result in the output
		//! \param[in] inputs input array for which we want to compute the average
		//! \param[out] output Array that will contain the average (in position 0)
		//! \param[in] stream CUDA stream (optional)
		void Average(VoxelType * output, cudaStream_t stream = nullptr) {
			DeviceArray<VoxelType> d_out(*temporaryBuffer);
			double multiplyFactor = 1.0 / temporaryBuffer->Length();
			Sum(temporaryBuffer->Pointer(), d_out.Pointer(), temporaryBuffer->Length(), (VoxelType) multiplyFactor, stream);
			checkCudaErrors(cudaMemcpy(output,d_out.Pointer(),sizeof(VoxelType),cudaMemcpyDeviceToHost));
		}

		//! Computes the minimum of an input array, placing the result in the output
		//! \param[in] inputs input array for which we want to compute the minimum
		//! \param[out] output Array that will contain the minimum (in position 0)
		//! \param[in] stream CUDA stream (optional)
		void Min(VoxelType * output, cudaStream_t stream = nullptr) {
			DeviceArray<VoxelType> d_out(*temporaryBuffer);
			Min(temporaryBuffer->Pointer(), d_out.Pointer(), temporaryBuffer->Length(), stream);
			checkCudaErrors(cudaMemcpy(output,d_out.Pointer(),sizeof(VoxelType),cudaMemcpyDeviceToHost));
		}

		// currently not available
		//! Computes the minimum of an input matrix, placing the result in the output
		//! \param[in] inputs input matrix for which we want to compute the minimum
		//! \param[out] output Array that will contain the minimum (in position 0)
		//! \param[in] stream CUDA stream (optional)
		//void Min(DeviceMatrix<VoxelType> & inputs, DeviceArray<VoxelType> & output, cudaStream_t stream = nullptr) {
		//	Min(inputs.Pointer(), output.Pointer(), inputs.Elements(), stream);
		//}

		//! Computes the minimum of an input array as well as its index within the array
		//! \param[in] inputs input array for which we want to compute the minimum
		//! \param[out] min Array that will contain the minimum (in position 0)
		//! \param[out] minIndex Array that will contain the index of the minimum within the array (in position 0)
		//! \param[in] stream CUDA stream (optional)
		void MinIndex(VoxelType * output, int * outidx, cudaStream_t stream = nullptr) {
			DeviceArray<int> minIndex( temporaryBuffer->Length());
			DeviceArray<VoxelType> d_out(*temporaryBuffer);
			MinIndex(temporaryBuffer->Pointer(), d_out.Pointer(), minIndex.Pointer(), temporaryBuffer->Length(), stream);
			checkCudaErrors(cudaMemcpy(output, d_out.Pointer(),sizeof(VoxelType),cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(outidx, minIndex.Pointer(),sizeof(int),cudaMemcpyDeviceToHost));
		}

		// currently not available
		//! Computes the minimum of an input matrix as well as its (1-D) index within the matrix
		//! \param[in] inputs input matrix for which we want to compute the minimum
		//! \param[out] min Array that will contain the minimum (in position 0)
		//! \param[out] minIndex Array that will contain the index of the minimum within the array (in position 0)
		//! \param[in] stream CUDA stream (optional)
		//void MinIndex(DeviceMatrix<VoxelType> & inputs, DeviceArray<VoxelType> & min, DeviceArray<int> & minIndex, cudaStream_t stream = nullptr) {
		//	MinIndex(inputs.Pointer(), min.Pointer(), minIndex.Pointer(), inputs.Elements(), stream);
		//}

		//! Computes the maximum of an input array, placing the result in the output
		//! \param[in] inputs input array for which we want to compute the maximum
		//! \param[out] output Array that will contain the maximum (in position 0)
		//! \param[in] stream CUDA stream (optional)
		void Max(VoxelType * output, cudaStream_t stream = nullptr) {
			DeviceArray<VoxelType> d_out(*temporaryBuffer);
 			Max(temporaryBuffer->Pointer(), d_out.Pointer(), temporaryBuffer->Length(), stream);
			checkCudaErrors(cudaMemcpy(output, d_out.Pointer(),sizeof(VoxelType),cudaMemcpyDeviceToHost));
 		}

		// currently not available
		//! Computes the maximum of an input matrix, placing the result in the output
		//! \param[in] inputs input matrix for which we want to compute the maximum
		//! \param[out] output Array that will contain the maximum (in position 0)
		//! \param[in] stream CUDA stream (optional)
		//void Max(DeviceMatrix<VoxelType> & inputs, DeviceArray<VoxelType> & output, cudaStream_t stream = nullptr) {
		//	Max(inputs.Pointer(), output.Pointer(), inputs.Elements(), stream);
		//}

		//! Computes the maximum of an input array as well as its index within the array
		//! \param[in] inputs input array for which we want to compute the minimum
		//! \param[out] max Array that will contain the minimum (in position 0)
		//! \param[out] maxIndex Array that will contain the index of the minimum within the array (in position 0)
		//! \param[in] stream CUDA stream (optional)
		void MaxIndex(VoxelType * output, int * outidx, cudaStream_t stream = nullptr) {
			DeviceArray<int> maxIndex( temporaryBuffer->Length());
			DeviceArray<VoxelType> d_out(*temporaryBuffer);
			MaxIndex(temporaryBuffer->Pointer(), d_out.Pointer(), maxIndex.Pointer(), temporaryBuffer->Length(), stream);
			checkCudaErrors(cudaMemcpy(output,d_out.Pointer(),sizeof(VoxelType),cudaMemcpyDeviceToHost));
			checkCudaErrors(cudaMemcpy(outidx,maxIndex.Pointer(),sizeof(int),cudaMemcpyDeviceToHost));
		}

		// currently not available
		//! Computes the maximum of an input matrix as well as its (1-D) index within the array
		//! \param[in] inputs input matrix for which we want to compute the minimum
		//! \param[out] max Array that will contain the minimum (in position 0)
		//! \param[out] maxIndex Array that will contain the index of the minimum within the array (in position 0)
		//! \param[in] stream CUDA stream (optional)
		//void MaxIndex(DeviceMatrix<VoxelType> & inputs, DeviceArray<VoxelType> & max, DeviceArray<int> & maxIndex, cudaStream_t stream = nullptr) {
		//	MaxIndex(inputs.Pointer(), max.Pointer(), maxIndex.Pointer(), inputs.Elements(), stream);
		//}
};

//! @}

#endif