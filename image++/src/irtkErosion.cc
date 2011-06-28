/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkErosion.h>

template <class VoxelType> irtkErosion<VoxelType>::irtkErosion()
{
	// Default connectivity.
	this->_Connectivity = CONNECTIVITY_26;
}

template <class VoxelType> irtkErosion<VoxelType>::~irtkErosion(void)
{}

template <class VoxelType> bool irtkErosion<VoxelType>::RequiresBuffering(void)
{
  return true;
}

template <class VoxelType> const char *irtkErosion<VoxelType>::NameOfClass()
{
  return "irtkErosion";
}

template <class VoxelType> void irtkErosion<VoxelType>::Initialize()
{
  // Do the initial set up
  this->irtkImageToImage<VoxelType>::Initialize();

  this->_offsets.Initialize(this->_input, this->_Connectivity);
}

template <class VoxelType> void irtkErosion<VoxelType>::Run()
{
  int i, x, y, z, t, maskSize;
  VoxelType value;
  VoxelType *ptr2current, *ptr2offset;

  // Do the initial set up
  this->Initialize();

  maskSize = this->_offsets.GetSize();

  for (t = 0; t < this->_input->GetT(); t++) {
    for (z = 0; z < this->_input->GetZ(); z++) {
      for (y = 0; y < this->_input->GetY(); y++) {
        for (x = 0; x < this->_input->GetX(); x++) {
          if ((x == 0) || (x == this->_input->GetX()-1) ||
              (y == 0) || (y == this->_input->GetY()-1) ||
              (z == 0) || (z == this->_input->GetZ()-1)) {
            this->_output->Put(x, y, z, t, this->_input->Get(x, y, z, t));
          } else {
            value = this->_input->Get(x, y, z, t);
          	ptr2current = this->_input->GetPointerToVoxels(x, y, z, t);
          	for (i = 0; i < maskSize; ++i){
          		ptr2offset = ptr2current + this->_offsets(i);
          		if (*ptr2offset < value)
          			value = *ptr2offset;
          	}
            this->_output->Put(x, y, z, t, value);
          }
        }
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();
}


template class irtkErosion<irtkBytePixel>;
template class irtkErosion<irtkGreyPixel>;
template class irtkErosion<irtkRealPixel>;
