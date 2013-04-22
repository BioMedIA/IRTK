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
	d_input = _input;
}

template <class VoxelType> 
void irtkCUMaskFilter<VoxelType>::setOutput(irtkCUGenericImage<VoxelType>* _output)
{
	d_output = _output;
}


template <class VoxelType> 
irtkCUGenericImage<VoxelType>* irtkCUMaskFilter<VoxelType>::getOuput()
{
	run();
	return d_output;
}

template <class VoxelType> 
void irtkCUMaskFilter<VoxelType>::run()
{
	if(d_output == NULL)
		d_output = new irtkCUGenericImage<VoxelType>(*d_input);

	mask_image_gpu(d_input->GetX()*d_input->GetY()*d_input->GetZ(), 
		d_input->getDevicePoiner(), d_output->getDevicePoiner(), d_mask->getDevicePoiner(), nullValue, maskValue);
}


template class irtkCULib_DLLAPI irtkCUMaskFilter<char>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<unsigned char>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<short>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<unsigned short>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<int>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<unsigned int>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<float>;
template class irtkCULib_DLLAPI irtkCUMaskFilter<double>;