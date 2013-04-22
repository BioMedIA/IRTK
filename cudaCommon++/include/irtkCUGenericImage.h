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

#pragma once

#include <iostream>
#include <stdlib.h>

#include <irtkObject.h>
#include <irtkImage.h>
#include <irtkBaseImage.h>
#include <irtkGenericImage.h>

#include "irtkCUSharedLibMacros.h"

#ifdef USE_CUDA_NPP
#include <nppdefs.h>
#endif

template <class T> class IRTKCUVolume3D;

template <class VoxelType> class irtkCULib_DLLAPI irtkCUGenericImage : public irtkGenericImage<VoxelType>
{
protected:

	/// Pointer to image data on GPU
	// TODO:
	// provide also fast read only texture access
	VoxelType* d_matrix;

public:
	/// Default constructor
	irtkCUGenericImage();

	/// Constructor from image file
	irtkCUGenericImage(char *);

	/// Constructor for given image size
	irtkCUGenericImage(int, int, int, int = 1);

	/// Copy constructor for image 
	irtkCUGenericImage(const irtkCUGenericImage &);

	/// Constructor for given image attributes
	irtkCUGenericImage(const irtkImageAttributes &);

	/// Copy constructor for image of different type
	template <class T> irtkCUGenericImage(const irtkGenericImage<T> &);

	/// Copy constructor for image of different type
	template <class T> irtkCUGenericImage(const irtkCUGenericImage<T> &);

	/// Destructor
	~irtkCUGenericImage(void);

	/// Get Pointer on GPU device
	VoxelType* getDevicePoiner() const;

	/// Function for pixel access via pointers
	VoxelType *GetPointerToVoxels(int = 0, int = 0, int = 0, int = 0) const;

	/// Function to convert pixel to index
	int VoxelToIndex(int, int, int, int = 0) const;

	/// Function for pixel get access
	VoxelType Get(int, int, int, int = 0) const;

	/// Function for pixel access from via operators
	VoxelType& operator()(int, int, int, int = 0);

	/// TODO other convenience functions...
	
	// GPU function
	/// Minimum and maximum pixel values get accessors using NPP
	void GetMinMax(VoxelType *min, VoxelType *max) const;  

	//Device function
	/// Minimum and maximum pixel values get accessors using NPP
	VoxelType getAverage(int = 1) const; 

#ifdef USE_CUDA_NPP
	//Device function
	/// Minimum and maximum pixel values get accessors using NPP
	//NPP only!
	void GetHistogram(int*  hist, int numBins) const; 
#endif

	//Device function
	// Normalizes image depending on type
	// needs image maximum and minimum values 
	// min and max can by computed with GetMinMax
	// if type == float || double: normalization between [-1,1]
	// else normalization to type range
	void normalize(VoxelType minVal, VoxelType maxVal) const; 

private:
	void cudaCheckSets();

	// internal performance primitive implementations of minmax
	// for different datatypes and the use of NPP
	void GetMinMax_internal(VoxelType *min, VoxelType *max) const;  

protected:
#ifdef USE_CUDA_NPP
	//internal histogram computation
	void computeHistogram_internal(int*  hist, int numBins, NppiSize oSizeROI, Npp32s *histDevice, 
		const int levelCount, VoxelType max_, VoxelType min_) const;
#endif
};


template <class VoxelType> inline VoxelType *irtkCUGenericImage<VoxelType>::GetPointerToVoxels(int x, int y, int z, int t) const
{
#ifdef __CUDACC__
	return &d_matrix[t*z*y*x];
#else
#ifdef NO_BOUNDS
  return &(this->_matrix[t][z][y][x]);
#else
  if ((x >= this->_attr._x) || (x < 0) || (y >= this->_attr._y) || (y < 0) || (z >= this->_attr._z) || (z < 0) || (t >= this->_attr._t) || (t < 0)) {
    cout << "irtkGenericImage<Type>::GetPointerToVoxels: parameter out of range\n";
    cout << x << " " << y << " " << z << " " << t << endl;
    return NULL;
  } else {
    return &(this->_matrix[t][z][y][x]);
  }
#endif
#endif
}


template <class VoxelType> inline int irtkCUGenericImage<VoxelType>::VoxelToIndex(int x, int y, int z, int t) const
{
#ifdef __CUDACC__
	return (&(d_matrix[t*z*y*x]) - &(d_matrix[0]));
#else
#ifdef NO_BOUNDS
	return (&(this->_matrix[t][z][y][x]) - &(this->_matrix[0][0][0][0]));
#else
	if ((x >= this->_attr._x) || (x < 0) || (y >= this->_attr._y) || (y < 0) || (z >= this->_attr._z) || (z < 0) || (t >= this->_attr._t) || (t < 0)) {
		cout << "irtkGenericImage<Type>::VoxelToIndex: parameter out of range\n";
		return 0;
	} else {
		return (&(this->_matrix[t][z][y][x]) - &(this->_matrix[0][0][0][0]));
	}
#endif
#endif
}


template <class VoxelType> inline VoxelType irtkCUGenericImage<VoxelType>::Get(int x, int y, int z, int t) const
{
#ifdef __CUDACC__
	return d_matrix[t*z*y*x];
#else
#ifdef NO_BOUNDS
	return (this->_matrix[t][z][y][x]);
#else
	if ((x >= this->_attr._x) || (x < 0) || (y >= this->_attr._y) || (y < 0) || (z >= this->_attr._z) || (z < 0) || (t >= this->_attr._t) || (t < 0)) {
		cout << "irtkGenericImage<Type>::Get: parameter out of range\n";
		return 0;
	} else {
		return(this->_matrix[t][z][y][x]);
	}
#endif
#endif
}

template <class VoxelType> inline VoxelType& irtkCUGenericImage<VoxelType>::operator()(int x, int y, int z, int t)
{
#ifdef __CUDACC__
	return d_matrix[t*z*y*x];
#else
#ifdef NO_BOUNDS
	return (this->_matrix[t][z][y][x]);
#else
	if ((x >= this->_attr._x) || (x < 0) || (y >= this->_attr._y) || (y < 0) || (z >= this->_attr._z) || (z < 0) || (t >= this->_attr._t) || (t < 0)) {
		cout << "irtkGenericImage<Type>::(): parameter out of range\n";
		return this->_matrix[0][0][0][0];
	} else {
		return (this->_matrix[t][z][y][x]);
	}
#endif
#endif
}
