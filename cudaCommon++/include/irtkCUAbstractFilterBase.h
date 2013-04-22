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

#ifndef _irtkCUFilterBase_H

#define _irtkCUFilterBase_H


#include <iostream>
#include <stdlib.h>

#include <irtkObject.h>
#include <irtkCUImage.h>
#include <irtkBaseImage.h>
#include <irtkGenericImage.h>
#include <irtkCUGenericImage.h>

#include "irtkCUSharedLibMacros.h"

#include <vector>
using namespace std;


/*

Abstract base class for filters using cuda

*/

template <class VoxelType> irtkCULib_DLLAPI  class irtkCUAbstractFilterBase : public irtkObject
{
protected: 
		irtkCUAbstractFilterBase() {};
		irtkCUGenericImage<VoxelType>* d_input; 
		irtkCUGenericImage<VoxelType>* d_output; 

		// convinience -- filters are usually evaluated when output is requested
		virtual void run() = 0;
public:

	// sets the input
	virtual void setInput(irtkCUGenericImage<VoxelType>* input) = 0;

	// sets the output -- if optional allocate pointer in getOutput
	virtual void setOutput(irtkCUGenericImage<VoxelType>* output) = 0;

	//executes the filter and returns the output
	virtual irtkCUGenericImage<VoxelType>* getOuput() = 0;

};


#endif