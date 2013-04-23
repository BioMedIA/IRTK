/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : 
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : 
  Version   : 
  Changes   : $Author: bkainz $

=========================================================================*/

#ifndef _irtkCUMaskFilter_H

#define _irtkCUMaskFilter_H


#include <iostream>
#include <stdlib.h>

#include "irtkCUAbstractFilterBase.h"

#include "irtkCUSharedLibMacros.h"
#include "irtkCUImage.h"

#include <vector>
using namespace std;


/*

Masking an image using cuda

*/

template <class VoxelType> irtkCULib_DLLAPI
class irtkCUMaskFilter : public irtkCUAbstractFilterBase<VoxelType>
{
private:
	irtkCUGenericImage<VoxelType>* d_mask; 
	VoxelType nullValue, maskValue;

public:
	irtkCUMaskFilter();
	~irtkCUMaskFilter();

	void setMask(irtkCUGenericImage<VoxelType>* _mask, VoxelType _nullValue = 0, VoxelType _maskValue = 0);
	void setInput(irtkCUGenericImage<VoxelType>* _input);
	void setOutput(irtkCUGenericImage<VoxelType>* _output);
	irtkCUGenericImage<VoxelType>* getOuput();
	void run();

};

/// Unsigned char image
typedef class irtkCUMaskFilter<irtkBytePixel> irtkByteFilter;
/// Short image
typedef class irtkCUMaskFilter<irtkGreyPixel> irtkGreyFilter;
/// Float image
typedef class irtkCUMaskFilter<irtkRealPixel> irtkRealFilter;

#endif
