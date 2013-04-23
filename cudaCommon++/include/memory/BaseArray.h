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

#ifndef GPUMLib_BaseArray_h
#define GPUMLib_BaseArray_h

#ifndef HAVE_CXX11_NULLPTR
#include "nullptr.h"
#endif

//! \addtogroup memframework Host (CPU) and device (GPU) memory access framework
//! @{

template <class Type> class CudaArray;

//! Base class for HostArray and DeviceArray classes (Array base class)
template <class Type> class BaseArray {
	friend class CudaArray<Type>;

	protected:
		Type * arrayData;
		int size;

		void Init() {
			arrayData = nullptr;
			size = 0;
		}

		BaseArray() {
			Init();
		}

		#ifdef Cx11
		void MoveFrom(BaseArray<Type> & other) {
			size = other.size;
			arrayData = other.arrayData;

			other.Init();
		}
		#endif

		virtual void Alloc(int size) = 0;

	public:
		//! Disposes the array.
		virtual void Dispose() = 0;

		//! Gets the length of the array. You can use this function to check if the array was effectively allocated.
		//! \return the number of elements of the array
		int Length() const {
			return size;
		}

		//! Gets a pointer to the array data
		//! \attention Use with caution
		//! \return a pointer to the array data
		Type * Pointer() const {
			return (size > 0) ? arrayData : nullptr;
		}

		//! Resizes the array without preserving its data
		//! \param size new size of the array
		//! \return the number of elements of the array after being resized.
		int ResizeWithoutPreservingData(int size) {
			if (size != this->size) {
				Dispose();
				Alloc(size);
			}

			return Length();
		}
};

//! @}

#endif