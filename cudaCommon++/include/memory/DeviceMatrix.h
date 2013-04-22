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


#ifndef GPUMLib_DeviceMatrix_h
#define GPUMLib_DeviceMatrix_h

#include <cublas.h>

#include "../reduction/CudaDefinitions.h"
#include "HostMatrix.h"

//! \addtogroup memframework Host (CPU) and device (GPU) memory access framework
//! @{

//! Create a matrix of any type, on the device, that automatically manages the memory used to hold its elements
template <class Type> class DeviceMatrix : public BaseMatrix<Type> {
	private:
		void Alloc(int rows, int columns) {
			assert(rows > 0 && columns > 0);

			int elements = rows * columns;

			if (cudaMalloc((void **) &(this->matrixData), elements * sizeof(Type)) == cudaSuccess) {
				this->rows = rows;
				this->columns = columns;
			} else {
				this->Init();
			}
		}

		void Assign(const HostMatrix<Type> & other) {
			this->storingOrder = (other.IsRowMajor()) ? RowMajor : ColumnMajor;
			int elements = this->ResizeWithoutPreservingData(other.Rows(), other.Columns());
			if (elements > 0) cudaMemcpy(this->matrixData, other.Pointer(), elements * sizeof(Type), cudaMemcpyHostToDevice);
		}

		void Assign(const DeviceMatrix<Type> & other) {
			this->storingOrder = other.storingOrder;
			int elements = this->ResizeWithoutPreservingData(other.rows, other.columns);
			if (elements > 0) cudaMemcpy(this->matrixData, other.Pointer(), elements * sizeof(Type), cudaMemcpyDeviceToDevice);
		}

	public:
		void Dispose() {
			if (this->matrixData != nullptr) cudaFree(this->matrixData);
			this->Init();
		}

		//! Constructs an empty matrix
		//! \param storingOrder defines if the matrix uses the row-major or column-major order to store the information
		DeviceMatrix(StoringOrder storingOrder = RowMajor) : BaseMatrix<Type>(storingOrder) {}

		//! Constructs a matrix with a given number of rows and columns
		//! \param rows the number of rows
		//! \param columns the number of columns
		//! \param storingOrder defines if the matrix uses the row-major or column-major order to store the information
		DeviceMatrix(int rows, int columns, StoringOrder storingOrder = RowMajor) : BaseMatrix<Type>(storingOrder) {
			assert(rows > 0 && columns > 0);
			this->ResizeWithoutPreservingData(rows, columns);
		}

		//! Constructs a matrix identical to another
		//! \param other another matrix
		DeviceMatrix(const DeviceMatrix<Type> & other) {
			Assign(other);
		}

		//! Constructs a matrix identical to an host matrix
		//! \param other host matrix
		DeviceMatrix(const HostMatrix<Type> & other) {
			Assign(other);
		}

		#ifdef Cx11
		//! Constructs a matrix using the elements of a device temporary matrix (rvalue)
		//! \param temporaryMatrix temporary device matrix containing the elements
		DeviceMatrix(DeviceMatrix<Type> && temporaryMatrix) {
			this->MoveFrom(temporaryMatrix);
		}
		#endif

		//! Destructor
		~DeviceMatrix() {
			Dispose();
		}

		//! Transforms this matrix into an matrix identical to an host matrix
		//! \param other host matrix	
		//! \return a reference to this matrix
		//! \attention The storing order (major-row or major-column) becomes the same of the other matrix.
		//! \sa IsRowMajor
		DeviceMatrix<Type> & operator = (const HostMatrix<Type> & other) {
			Assign(other);
			return *this;
		}

		//! Transforms this matrix into an matrix identical to the other
		//! \param other other matrix
		//! \return a reference to this matrix
		//! \attention The storing order (major-row or major-column) becomes the same of the other matrix.
		//! \sa IsRowMajor
		DeviceMatrix<Type> & operator = (const DeviceMatrix<Type> & other) {
			Assign(other);
			return *this;
		}

		#ifdef Cx11
		//! Replaces this matrix using a temporary matrix (rvalue)
		//! \param temporaryMatrix temporary matrix
		//! \return a reference to this matrix
		DeviceMatrix<Type> & operator = (const DeviceMatrix<Type> && temporaryMatrix) {
			if (this != &temporaryMatrix) {
				Dispose();
				this->MoveFrom(temporaryMatrix);
			}
			
			return *this;
		}
		#endif

		//! Releases its own resources (elements) and obtains ownership of another matrix resources. 
		//! The other matrix will no longer have any elements. 
		//! In other words, it moves the elements from one device matrix to another.
		//! \param other matrix containing the elements to be moved.
		void TransferOwnerShipFrom(DeviceMatrix<Type> & other) {
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
		DeviceMatrix<Type> Transpose() {
			HostMatrix<Type> transpose(*this);
			transpose.ReplaceByTranspose();

			return transpose;
		}

		//! Multiplies the matrix by its own transpose and places the result in matrix C. More specifically C = alpha * (A * transpose(A)) + beta * C. This method uses the CUBLAS library.
		//! \attention Matrix C is returned in column-major order. If beta is different than zero, Matrix C must already be in column-major order.
		//! \param C matrix C. The number of rows and columns of C must be identical to the number of rows of A.
		//! \param alpha alpha scalar parameter (default value is 1.0). 
		//! \param beta beta scalar parameter (default value is 0.0). If beta is different than zero, Matrix C must already be in column-major order.
		//! \sa IsRowMajor
		//! \warning The CUBLAS library must be initialized prior to the use of this method. Use cublasInit() to initialize the CUBLAS library. Don't forget to call cublasShutdown() when you no longer need to use the CUBLAS library.
		template <class VoxelType>
		void MultiplyBySelfTranspose(DeviceMatrix<Type> & C, VoxelType alpha = 1.0, VoxelType beta = 0.0) {
			assert(C.rows == this->rows && C.columns == this->rows);
			
			if (C.IsRowMajor()) {
				assert(beta == 0.0);
				C.storingOrder = ColumnMajor;
			}

			int ldAB = this->IsRowMajor() ? this->columns : this->rows;
			cublasSgemm(this->IsRowMajor() ? 'T' : 'N', this->IsRowMajor() ? 'N' : 'T', C.rows, C.columns, this->columns, alpha, this->matrixData, ldAB, this->matrixData, ldAB, beta, C.matrixData, C.rows);
		}

		//! Multiplies matrix A by Matrix B and places the result in C. More specifically C = alpha * (A * B) + beta * C. This method uses the CUBLAS library.
		//! \attention Matrix C is returned in column-major order. If beta is different than zero, Matrix C must already be in column-major order.
		//! \attention Best performance is achieved when all matrices are in column-major order.
		//! \param A matrix A. The number of columns of A must be identical to the number of B rows.
		//! \param B matrix B. The number of rows of B must be identical to the number of A columns.
		//! \param C matrix C. The number of rows of C must be identical to the number of rows of A and the number of columns of C must be identical to the number of columns of B.
		//! \param alpha alpha scalar parameter (default value is 1.0). 
		//! \param beta beta scalar parameter (default value is 0.0). If beta is different than zero, Matrix C must already be in column-major order.
		//! \sa IsRowMajor, ReplaceByTranspose, Transpose
		//! \warning The CUBLAS library must be initialized prior to the use of this method. Use cublasInit() to initialize the CUBLAS library. Don't forget to call cublasShutdown() when you no longer need to use the CUBLAS library.
		template <class VoxelType>
		static void Multiply(DeviceMatrix<Type> & A, DeviceMatrix<Type> & B, DeviceMatrix<Type> & C, VoxelType alpha = 1.0, VoxelType beta = 0.0) {
			assert(A.columns == B.rows && C.rows == A.rows && C.columns == B.columns);
			
			if (C.IsRowMajor()) {
				assert(beta == 0.0);
				C.storingOrder = ColumnMajor;
			}

			cublasSgemm(A.IsRowMajor() ? 'T' : 'N', B.IsRowMajor() ? 'T' : 'N', C.rows, C.columns, A.columns, alpha, A.matrixData, A.IsRowMajor() ? A.columns : A.rows, B.matrixData, B.IsRowMajor() ? B.columns : B.rows, beta, C.matrixData, C.rows);
		}
};

//! @}

#endif
