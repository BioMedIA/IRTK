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

#ifndef GPUMLib_BaseMatrix_h
#define GPUMLib_BaseMatrix_h

//! \addtogroup memframework Host (CPU) and device (GPU) memory access framework
//! @{

//! Defines the methods for storing information in the matrix (either row-major or column-major).
typedef enum {
	RowMajor,
	ColumnMajor
} StoringOrder;

template <class Type> class CudaMatrix;

//! Base class for HostMatrix and DeviceMatrix classes (Matrix base class)
template <class Type> class BaseMatrix {
	friend class CudaMatrix<Type>;

	protected:
		Type * matrixData;
		int rows;
		int columns;
		StoringOrder storingOrder;

		void Init() {
			matrixData = nullptr;
			rows = columns = 0;
		}

		BaseMatrix(StoringOrder storingOrder = RowMajor) {
			this->storingOrder = storingOrder;
			Init();
		}

		#ifdef Cx11
		void MoveFrom(BaseMatrix<Type> & other) {
			matrixData = other.matrixData;
			rows = other.rows;
			columns = other.columns;
			storingOrder = other.storingOrder;

			other.Init();
		}
		#endif

		virtual void Alloc(int rows, int columns) = 0;

	public:
		//! Disposes the matrix
		virtual void Dispose() = 0;

		//! Gets the number of rows of the matrix
		//! \return the number of rows of the matrix
		int Rows() const {
			return rows;
		}

		//! Gets the number of columns of the matrix
		//! \return the number of columns of the matrix
		int Columns() const {
			return columns;
		}

		//! Gets a pointer to the matrix data
		//! \attention Use with caution. Special attention should be given to how the matrix information is stored (row-major or column-major).
		//! \return a pointer to the matrix data
		//! \sa IsRowMajor
		Type * Pointer() const {
			return (Elements() > 0) ? matrixData : nullptr;
		}

		//! Gets the number of elements contained in the matrix
		//! \return the number of elements contained in the matrix
		int Elements() const {
			return rows * columns;
		}

		//! Indicates if the information in the matrix is stored in row-major order.
		//! \return True if the matrix information is stored in row-major order. False if the information is stored in column-major format.
		bool IsRowMajor() const {
			return (storingOrder == RowMajor);
		}

		//! Replaces this matrix by its transpose
		//! \attention This method is very fast, however it changes the method for storing information in the matrix (row-major or column-major).
		//! \sa IsRowMajor
		void ReplaceByTranspose() {
			int newRows = columns;

			columns = rows;
			rows = newRows;
			storingOrder = (IsRowMajor()) ? ColumnMajor : RowMajor;
		}

		//! Resizes the matrix without preserving its data
		//! \param rows the new number of rows
		//! \param columns the new number of columns
		//! \return the number of elements of the matrix after being resized.
		int ResizeWithoutPreservingData(int rows, int columns) { 
			int newElements = rows * columns;

			if (newElements != Elements()) {
				Dispose();
				Alloc(rows, columns);
			} else {
				this->rows = rows;
				this->columns = columns;
			}

			return Elements();
		}
};

//! @}

#endif