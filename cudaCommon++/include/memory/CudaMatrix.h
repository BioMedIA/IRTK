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

#ifndef GPUMLib_CudaMatrix_h
#define GPUMLib_CudaMatrix_h

#include "DeviceMatrix.h"

//! \addtogroup memframework Host (CPU) and device (GPU) memory access framework
//! @{

//! Create a matrix of any type, that automatically manages the memory used to hold its elements (data will be stored both on the host and on the device). \attention The data on the host might differ from the data on the device. Use UpdateDevice and UpdateHost to synchronize data.
template <class Type> class CudaMatrix {
	private:
		HostMatrix<Type> h;
		DeviceMatrix<Type> d;

		#ifdef Cx11
		void MoveFrom(BaseMatrix<Type> & other) {
			h.MoveFrom(other.h);
			d.MoveFrom(other.d);
		}
		#endif

	public:
		//! Constructs an empty matrix
		//! \param storingOrder defines if the matrix uses the row-major or column-major order to store the information
		CudaMatrix(StoringOrder storingOrder = RowMajor) : h(storingOrder), d(storingOrder) { }

		//! Constructs a matrix with a given number of rows and columns
		//! \param rows the number of rows
		//! \param columns the number of columns
		//! \param storingOrder defines if the matrix uses the row-major or column-major order to store the information
		CudaMatrix(int rows, int columns, StoringOrder storingOrder = RowMajor) : h(rows, columns, storingOrder), d(rows, columns, storingOrder) { }

		//! Constructs a matrix identical to another matrix
		//! \param other another matrix
		CudaMatrix(const CudaMatrix<Type> & other) : h(other.h), d(other.d) { }

		//! Constructs a matrix identical to an host matrix
		//! \param other another matrix
		CudaMatrix(const HostMatrix<Type> & other) : h(other), d(other) { }

		//! Constructs a matrix identical to a device matrix
		//! \param other device matrix
		CudaMatrix(const DeviceMatrix<Type> & other) : h(other), d(other) { }

		#ifdef Cx11
		//! Constructs a matrix using the elements of a temporary matrix (rvalue)
		//! \param temporaryMatrix temporary matrix containing the elements
		CudaMatrix(CudaMatrix<Type> && temporaryMatrix) {
			MoveFrom(temporaryMatrix);
		}
		#endif

		//! Transforms this matrix into an matrix identical to a host matrix
		//! \param other host matrix	
		//! \return a reference to this matrix
		//! \attention The storing order (major-row or major-column) becomes the same of the other matrix.
		//! \sa IsRowMajor
		CudaMatrix<Type> & operator = (const HostMatrix<Type> & other) {
			h = other;
			d = other;

			return *this;
		}

		//! Transforms this matrix into a matrix identical to a device matrix
		//! \param other device matrix
		//! \return a reference to this matrix
		//! \attention The storing order (major-row or major-column) becomes the same of the other matrix.
		//! \sa IsRowMajor
		CudaMatrix<Type> & operator = (const DeviceMatrix<Type> & other) {
			h = other;
			d = other;

			return *this;
		}

		//! Transforms this matrix into a matrix identical to another
		//! \param other other matrix
		//! \return a reference to this matrix
		//! \attention The storing order (major-row or major-column) becomes the same of the other matrix.
		//! \sa IsRowMajor
		CudaMatrix<Type> & operator = (const CudaMatrix<Type> & other) {
			h = other.h;
			d = other.d;

			return *this;
		}

		#ifdef Cx11
		//! Replaces this matrix using a temporary matrix (rvalue)
		//! \param temporaryMatrix temporary matrix
		//! \return a reference to this matrix
		//! \attention The storing order (major-row or major-column) becomes the same of the other matrix.
		//! \sa IsRowMajor
		CudaMatrix<Type> & operator = (const CudaMatrix<Type> && temporaryMatrix) {
			if (this != &temporaryMatrix) {
				Dispose();
				MoveFrom(temporaryMatrix);
			}
			
			return *this;
		}
		#endif

		//! Updates the device matrix data with the host matrix data
		//! \attention The storing order (major-row or major-column) of the device matrix becomes the same as the host matrix.
		//! \sa IsRowMajor
		void UpdateDevice() {
			d = h;
		}

		//! Updates the host matrix data with the device matrix data
		//! \attention The storing order (major-row or major-column) of the host matrix becomes the same as the device matrix.
		//! \sa IsRowMajor
		void UpdateHost() {
			h = d;
		}

		//! Gets a pointer to the host matrix data
		//! \attention Use with caution. Special attention should be given to how the matrix information is stored (row-major or column-major).
		//! \return a pointer to the matrix data
		//! \sa IsRowMajor
		Type * HostPointer() const {
			return h.Pointer();
		}

		//! Gets a pointer to the device matrix data
		//! \attention Use with caution
		//! \return a pointer to the device matrix data
		Type * DevicePointer() const {
			return d.Pointer();
		}

		//! Gets a the device matrix
		//! \attention Use with caution
		//! \return The device matrix
		DeviceMatrix<Type> & GetDeviceMatrix() {
			return d;
		}

		//! Gets a the host matrix
		//! \attention Use with caution
		//! \return The host matrix
		HostMatrix<Type> & GetHostMatrix() {
			return h;
		}

		//! Gets a reference to an element of the host matrix
		//! \param row row of the desired element
		//! \param column column of the desired element
		//! \return a reference to an element desired, based on the row and column specified
		Type & operator()(int row, int column) {
			return h(row, column);
		}

		//! Gets an element of the host matrix
		//! \param row row of the desired element
		//! \param column column of the desired element
		//! \return the element desired, based on the row and column specified
		Type operator()(int row, int column) const {
			return h(row, column);
		}

		//! Gets the number of rows of the matrix
		//! \return the number of rows of the matrix
		int Rows() const {
			return h.Rows();
		}

		//! Gets the number of columns of the matrix
		//! \return the number of columns of the matrix
		int Columns() const {
			return h.Columns();
		}

		//! Gets the number of elements contained in the matrix
		//! \return the number of elements contained in the matrix
		int Elements() const {
			return h.Elements();
		}

		//! Indicates if the information in the matrix is stored in row-major order.
		//! \return True if the matrix information is stored in row-major order. False if the information is stored in column-major format.
		bool IsRowMajor() const {
			return h.IsRowMajor();
		}

		//! Replaces this matrix by its transpose
		//! \attention This method is very fast, however it changes the method for storing information in the matrix (row-major or column-major).
		//! \sa IsRowMajor
		void ReplaceByTranspose() {
			h.ReplaceByTranspose();
			d.ReplaceByTranspose();
		}

		//! Gets the transposed of the matrix
		//! \return the transposed of the matrix
		//! \attention The returned matrix does not use the same method (row-major or column-major) for storing information as this matrix.
		//! \sa ReplaceByTranspose, IsRowMajor
		CudaMatrix<Type> Transpose() {
			CudaMatrix<Type> transpose(h.Transpose());

			return transpose;
		}

		//! Resizes the matrix without preserving its data
		//! \param rows the new number of rows
		//! \param columns the new number of columns
		//! \return the number of elements of the matrix after being resized.
		int ResizeWithoutPreservingData(int rows, int columns) {
			int he = h.ResizeWithoutPreservingData(rows, columns);
			int de = d.ResizeWithoutPreservingData(rows, columns);

			return (he > de) ? de : he;
		}

		//! Disposes the matrix
		void Dispose() {			
			d.Dispose();
			h.Dispose();
		}

		//! //! Updates the host information and disposes the device matrix
		void DisposeDevice() {
			UpdateHost();
			d.Dispose();
		}

		//! Updates the device information and disposes the host matrix
		void DisposeHost() {
			UpdateDevice();
			h.Dispose();
		}
};

//! @}


#endif