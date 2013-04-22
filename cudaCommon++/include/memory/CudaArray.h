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

#ifndef GPUMLib_CudaArray_h
#define GPUMLib_CudaArray_h

#include "DeviceArray.h"

//! \addtogroup memframework Host (CPU) and device (GPU) memory access framework
//! @{

//! Create an array of any type, that automatically manages the memory used to hold its elements (data will be stored both on the host and on the device). \attention The data on the host might differ from the data on the device. Use UpdateDevice and UpdateHost to syncronize data.
template <class Type> class CudaArray {
	private:
		HostArray<Type> h;
		DeviceArray<Type> d;

		#ifdef Cx11
		void MoveFrom(BaseArray<Type> & other) {
			h.MoveFrom(other.h);
			d.MoveFrom(other.d);
		}
		#endif

	public:
		//! Constructs an array with no elements
		CudaArray() { }

		//! Constructs an array with size elements
		//! \param size number of elements of the array
		explicit CudaArray(int size) : h(size), d(size) { }

		//! Constructs an array with the same elements as another array
		//! \param originalArray original device array from where to copy the elements
		CudaArray(const DeviceArray<Type> & originalArray) : h(originalArray), d(originalArray) { }

		//! Constructs an array with the same elements as another array
		//! \param originalArray original array from where to copy the elements
		CudaArray(const HostArray<Type> & originalArray) : h(originalArray), d(originalArray) {	}

		#ifdef Cx11
		//! Constructs an array using the elements of a temporary array (rvalue)
		//! \param temporaryArray temporary array containing the elements
		CudaArray(CudaArray<Type> && temporaryArray) {
			MoveFrom(temporaryArray);
		}
		#endif

		//! Transforms this array into an array identical to another array
		//! \param originalArray original device array from where to copy the elements
		//! \return a reference to this array
		CudaArray<Type> & operator = (const CudaArray<Type> & originalArray) {
			h = originalArray.h;
			d = originalArray.d;
			
			return *this;
		}

		//! Transforms this array into an array identical to another array
		//! \param originalArray original device array from where to copy the elements
		//! \return a reference to this array
		CudaArray<Type> & operator = (const DeviceArray<Type> & originalArray) {
			d = originalArray;
			h = d;
			
			return *this;
		}

		//! Transforms this array into an array identical to another array
		//! \param originalArray original array from where to copy the elements
		//! \return a reference to this array
		CudaArray<Type> & operator = (const HostArray<Type> & originalArray) {
			h = originalArray;
			d = h;
			
			return *this;
		}

		#ifdef Cx11
		//! Replaces the elements of this array by the elements of a temporary array (rvalue)
		//! \param temporaryArray temporary array containing the elements
		//! \return a reference to this array
		CudaArray<Type> & operator = (const CudaArray<Type> && temporaryArray) {
			if (this != &temporaryArray) {
				Dispose();
				MoveFrom(temporaryArray);
			}
			
			return *this;
		}
		#endif

		//! Gets a reference to an element of the host array
		//! \param element position of the desired element
		//! \return a reference to an element desired
		Type & operator [] (int element) {
			return h[element];
		}

		//! Gets an element of the host array
		//! \param element position of the desired element
		//! \return the element desired
		Type operator [] (int element) const {
			return h[element];
		}

		//! Gets the lenght of the array
		//! \return the number of elements of the array
		int Length() const {
			return d.Length();
		}

		//! Gets a pointer to the host array data
		//! \attention Use with caution
		//! \return a pointer to the host array data
		Type * HostPointer() const {
			return h.Pointer();
		}

		//! Gets a pointer to the array data
		//! \attention Use with caution
		//! \return a pointer to the array data
		Type * DevicePointer() const {
			return d.Pointer();
		}

		//! Resizes the array without preserving its data
		//! \param size new size of the array
		//! \return the number of elements of the array after being resized.
		int ResizeWithoutPreservingData(int size) {
			int he = h.ResizeWithoutPreservingData(size);
			int de = d.ResizeWithoutPreservingData(size);

			return (he > de) ? de : he;
		}

		//! Updates the device array data with the host array data
		void UpdateDevice() {
			d = h;
		}

		//! Updates the host array data with the device array data
		void UpdateHost() {
			h = d;
		}

		//! Gets the device array
		//! \attention Use with caution
		//! \return The device array
		DeviceArray<Type> & GetDeviceArray() {
			return d;
		}

		//! Gets the device array
		//! \attention Use with caution
		//! \return The device array
		HostArray<Type> & GetHostArray() {
			return h;
		}

		//! Disposes the array
		void Dispose() {			
			d.Dispose();
			h.Dispose();
		}

		//! Updates the host information and disposes the device array
		void DisposeDevice() {
			UpdateHost();
			d.Dispose();
		}

		//! Updates the device information and disposes the Host array
		void DisposeHost() {
			UpdateDevice();
			h.Dispose();
		}
};

//! @}

#endif