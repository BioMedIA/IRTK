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

#ifndef GPUMLib_HostArray_h
#define GPUMLib_HostArray_h

#include <new>
#include <assert.h>
#include <cuda_runtime.h>

#include "BaseArray.h"

using namespace std;

template <class Type> class DeviceArray;

//! \addtogroup memframework Host (CPU) and device (GPU) memory access framework
//! @{

//! Create an array of any type, on the host, that automatically manages the memory used to hold its elements
template <class Type> class HostArray : public BaseArray<Type> {
	private:
		void Alloc(int size) {
			assert(size >= 0);
			
			if (size > 0) {
				this->arrayData = new (nothrow) Type[size];
				this->size = (this->arrayData != nullptr) ? size : 0;
			} else {
				this->Init();
			}
		}
		
	public:
		void Dispose() {
			if (this->size > 0) delete [] this->arrayData;
			this->Init();
		}

		//! Constructs an array with no elements
		HostArray() {}

		//! Constructs an array with size elements
		//! \param size number of elements of the array
		explicit HostArray(int size) {
			Alloc(size);
		}

		//! Constructs an array with the same elements as another array
		//! \param originalArray original device array from where to copy the elements
		HostArray(const DeviceArray<Type> & originalArray) {
			Alloc(originalArray.Length());
			if (this->size > 0) cudaMemcpy(this->arrayData, originalArray.Pointer(), this->size * sizeof(Type), cudaMemcpyDeviceToHost);
		}

		//! Constructs an array with the same elements as another array
		//! \param originalArray original array from where to copy the elements
		HostArray(const HostArray<Type> & originalArray) {
			Alloc(originalArray.Length());
			if (this->size > 0) memcpy(this->arrayData, originalArray.Pointer(), this->size * sizeof(Type));
		}

		#ifdef Cx11
		//! Constructs an array using the elements of a temporary array (rvalue)
		//! \param temporaryArray temporary array containing the elements
		HostArray(HostArray<Type> && temporaryArray) {
			this->MoveFrom(temporaryArray);
		}
		#endif

		//! Destructor
		~HostArray() {
			Dispose();
		}
		
		//! Transforms this array into an array identical to another array
		//! \param originalArray original device array from where to copy the elements
		//! \return a reference to this array
		HostArray<Type> & operator = (const DeviceArray<Type> & originalArray) {
			int newSize = originalArray.Length();
		
			this->ResizeWithoutPreservingData(newSize);
			if (this->size > 0) cudaMemcpy(this->arrayData, originalArray.Pointer(), this->size * sizeof(Type), cudaMemcpyDeviceToHost);
			
			return *this;
		}

		//! Transforms this array into an array identical to another array
		//! \param originalArray original array from where to copy the elements
		//! \return a reference to this array
		HostArray<Type> & operator = (const HostArray<Type> & originalArray) {
			int newSize = originalArray.Length();
		
			this->ResizeWithoutPreservingData(newSize);
			if (this->size > 0) memcpy(this->arrayData, originalArray.Pointer(), this->size * sizeof(Type));
			
			return *this;
		}

		#ifdef Cx11
		//! Replaces the elements of this array by the elements of a temporary array (rvalue)
		//! \param temporaryArray temporary array containing the elements
		//! \return a reference to this array
		HostArray<Type> & operator = (const HostArray<Type> && temporaryArray) {
			if (this != &temporaryArray) {
				Dispose();
				this->MoveFrom(temporaryArray);
			}
			
			return *this;
		}
		#endif

		//! Releases its own resources (elements) and obtains ownership of another array resources. 
		//! The other array will no longer have any elements. 
		//! In other words, it moves the elements from one device array to another.
		//! \param other array containing the elements to be moved.
		void TransferOwnerShipFrom(HostArray<Type> & other) {
			if (this != &other) {
				Dispose();
				this->size = other.size;
				this->arrayData = other.arrayData;
				other.Init();
			}
		}

		//! Gets a reference to an element of the array
		//! \param element position of the desired element
		//! \return a reference to an element desired
		Type & operator [] (int element) {
			assert(element < this->size);
			return this->arrayData[element];
		}

		//! Gets an element of the array
		//! \param element position of the desired element
		//! \return the element desired
		Type operator [] (int element) const {
			assert(element < this->size);
			return this->arrayData[element];
		}
};

//! @}
#endif
