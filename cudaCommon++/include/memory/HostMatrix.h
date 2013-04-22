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

#ifndef GPUMLib_HostMatrix_h
#define GPUMLib_HostMatrix_h

#include <new>
#include <assert.h>
#include <cuda_runtime.h>
#include <string.h>

#include "BaseMatrix.h"

using namespace std;

//! \addtogroup memframework Host (CPU) and device (GPU) memory access framework
//! @{

template <class Type> class DeviceMatrix;

//! Create a matrix of any type, on the host, that automatically manages the memory used to hold its elements
template <class Type> class HostMatrix : public BaseMatrix<Type> {
	private :
		void Alloc(int rows, int columns) {
			assert(rows > 0 && columns > 0);

			int elements = rows * columns;
			this->matrixData = new (nothrow) Type [elements];
			if (this->matrixData != nullptr) {
				this->rows = rows;
				this->columns = columns;
			} else {
				this->Init();
			}
		}

		void Assign(const HostMatrix<Type> & other) {
			this->storingOrder = other.storingOrder;
			int elements = this->ResizeWithoutPreservingData(other.rows, other.columns);
			if (elements > 0) memcpy(this->matrixData, other.matrixData, elements * sizeof(Type));
		}

		void Assign(const DeviceMatrix<Type> & other) {
			this->storingOrder = (other.IsRowMajor()) ? RowMajor : ColumnMajor;
			int elements = this->ResizeWithoutPreservingData(other.Rows(), other.Columns());
			if (elements > 0) cudaMemcpy(this->matrixData, other.Pointer(), elements * sizeof(Type), cudaMemcpyDeviceToHost);
		}

	public :
		void Dispose() {
			if (this->matrixData != nullptr) delete [] this->matrixData;
			this->Init();
		}

		//! Constructs an empty matrix
		//! \param storingOrder defines if the matrix uses the row-major or column-major order to store the information
		HostMatrix(StoringOrder storingOrder = RowMajor) : BaseMatrix<Type>(storingOrder) {}

		//! Constructs a matrix with a given number of rows and columns
		//! \param rows the number of rows
		//! \param columns the number of columns
		//! \param storingOrder defines if the matrix uses the row-major or column-major order to store the information
		HostMatrix(int rows, int columns, StoringOrder storingOrder = RowMajor) : BaseMatrix<Type>(storingOrder) {
			assert(rows > 0 && columns > 0);
			this->ResizeWithoutPreservingData(rows, columns);
		}

		//! Constructs a matrix identical to the other
		//! \param other another matrix
		HostMatrix(const HostMatrix<Type> & other) {
			Assign(other);
		}

		//! Constructs a matrix identical to a device matrix
		//! \param other device matrix
		HostMatrix(const DeviceMatrix<Type> & other) {
			Assign(other);
		}

		#ifdef Cx11
		//! Constructs a matrix using the elements of a temporary matrix (rvalue)
		//! \param temporaryMatrix temporary matrix containing the elements
		HostMatrix(HostMatrix<Type> && temporaryMatrix) {
			this->MoveFrom(temporaryMatrix);
		}
		#endif

		//! Destructor
		~HostMatrix() {
			Dispose();
		}

		//! Transforms this matrix into an matrix identical to the other
		//! \param other other matrix
		//! \return a reference to this matrix
		//! \attention The storing order (major-row or major-column) becames the same of the other matrix.
		//! \sa IsRowMajor
		HostMatrix<Type> & operator = (const HostMatrix<Type> & other) {
			Assign(other);
			return *this;
		}

		//! Transforms this matrix into an matrix identical a device matrix
		//! \param other device matrix
		//! \return a reference to this matrix
		//! \attention The storing order (major-row or major-column) becames the same of the other matrix.
		//! \sa IsRowMajor
		HostMatrix<Type> & operator = (const DeviceMatrix<Type> & other) {
			Assign(other);
			return *this;
		}

		#ifdef Cx11
		//! Replaces this matrix using a temporary matrix (rvalue)
		//! \param temporaryMatrix temporary matrix
		//! \return a reference to this matrix
		HostMatrix<Type> & operator = (const HostMatrix<Type> && temporaryMatrix) {
			if (this != &temporaryMatrix) {
				Dispose();
				this->MoveFrom(temporaryMatrix);
			}
			
			return *this;
		}
		#endif

		//! Releases its own resources (elements) and obtains ownership of another matrix resources. 
		//! The other matrix will no longer have any elements. 
		//! In other words, it moves the elements from one matrix to another.
		//! \param other matrix containing the elements to be moved.
		void TransferOwnerShipFrom(HostMatrix<Type> & other) {
			if (this != &other) {
				Dispose();
				this->matrixData = other.matrixData;
				this->rows = other.rows;
				this->columns = other.columns;
				this->storingOrder = other.storingOrder;

				other.Init();
			}
		}

		//! Gets the transposed of the matrix
		//! \return the transposed of the matrix
		//! \attention The returned matrix does not use the same method (row-major or column-major) for storing information as this matrix.
		//! \sa ReplaceByTranspose, IsRowMajor
		HostMatrix<Type> Transpose() {
			HostMatrix<Type> transpose(*this);
			transpose.ReplaceByTranspose();

			return transpose;
		}

		//! Gets a reference to an element of the matrix
		//! \param row row of the desired element
		//! \param column column of the desired element
		//! \return a reference to an element desired, based on the row and column specified
		Type & operator()(int row, int column) {
			assert(row >= 0 && row < this->rows && column >= 0 && column < this->columns);

			int index = (this->IsRowMajor()) ? row * this->columns + column : column * this->rows + row;

			return this->matrixData[index];
		}

		//! Gets an element of the matrix
		//! \param row row of the desired element
		//! \param column column of the desired element
		//! \return the element desired, based on the row and column specified
		Type operator()(int row, int column) const { 
			assert(row >= 0 && row < this->rows && column >= 0 && column < this->columns);

			int index = (this->IsRowMajor()) ? row * this->columns + column : column * this->rows + row;

			return this->matrixData[index];
		}
};

//! @}


#endif
