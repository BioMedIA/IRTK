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


#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

#include "../IO/irtkCUImage_kernels.cuh"
#include "irtkCUMaskFilter.h"

template <class VoxelType> 
irtkCUMaskFilter<VoxelType>::irtkCUMaskFilter()
{

}

template <class VoxelType> 
irtkCUMaskFilter<VoxelType>::~irtkCUMaskFilter()
{

}

template <class VoxelType> 
void irtkCUMaskFilter<VoxelType>::setMask(irtkCUGenericImage<VoxelType>* _mask, VoxelType _nullValue, VoxelType _maskValue)
{
	maskValue = _maskValue;
	nullValue = _nullValue;
	d_mask = _mask;
}

template <class VoxelType> 
void irtkCUMaskFilter<VoxelType>::setInput(irtkCUGenericImage<VoxelType>* _input)
{
	this->d_input = _input;
}

template <class VoxelType> 
void irtkCUMaskFilter<VoxelType>::setOutput(irtkCUGenericImage<VoxelType>* _output)
{
	this->d_output = _output;
}


template <class VoxelType> 
irtkCUGenericImage<VoxelType>* irtkCUMaskFilter<VoxelType>::getOuput()
{
	run();
	return this->d_output;
}

template <class VoxelType> 
void irtkCUMaskFilter<VoxelType>::run()
{
	if(this->d_output == NULL)
		this->d_output = new irtkCUGenericImage<VoxelType>(*this->d_input);

	mask_image_gpu(this->d_input->GetX()*this->d_input->GetY()*this->d_input->GetZ(), 
		this->d_input->getDevicePoiner(), this->d_output->getDevicePoiner(), d_mask->getDevicePoiner(), nullValue, maskValue);
}


template class irtkCULib_DLLAPI irtkCUMaskFilter<char>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<unsigned char>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<short>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<unsigned short>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<int>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<unsigned int>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<float>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<double>;
