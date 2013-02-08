/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2013-01-14$
  Changes   : $Author: cl6311 $

=========================================================================*/

#include <irtkImage.h>
#include <irtkMedianFilter.h>
#include <vector>
#include <algorithm>
#include <set>

template <class VoxelType> irtkMedianFilter<VoxelType>::irtkMedianFilter()
{
	// Default kernel radius.
	this->_kernelRadius = 5;

	/// irrelevant padding
	this->_mask = NULL;
}

template <class VoxelType> irtkMedianFilter<VoxelType>::~irtkMedianFilter(void)
{
}

template <class VoxelType> const char *irtkMedianFilter<VoxelType>::NameOfClass()
{
  return "irtkMedianFilter";
}

template <class VoxelType> void irtkMedianFilter<VoxelType>::Initialize()
{
  // Do the initial set up
  this->irtkImageToImage<VoxelType>::Initialize();
}

template <class VoxelType> bool irtkMedianFilter<VoxelType>::RequiresBuffering(void)
{
  return true;
}

template <class VoxelType> void irtkMedianFilter<VoxelType>::SetMask(irtkRealImage *mask)
{
  if (mask != NULL) {
    _mask = mask;
  } else {
    cerr << "irtkImageToImage::SetMask: mask is not an image\n";
    exit(1);
  }
}

template <class VoxelType> void irtkMedianFilter<VoxelType>::Run()
{
  int x, y, z, t;

  //int dim = 2*_kernelRadius + 1;

  // Do the initial set up
  this->Initialize();

  multiset<VoxelType> voxels;

  for (t = 0; t < this->_input->GetT(); t++) {
    for (z = 0; z < this->_input->GetZ(); z++) {
      for (y = 0; y < this->_input->GetY(); y++) {
        for (x = 0; x < this->_input->GetX(); x++) {
          if ((x < _kernelRadius) || (x > this->_input->GetX() - 1 - _kernelRadius) ||
              (y < _kernelRadius) || (y > this->_input->GetY() - 1 - _kernelRadius) ||
              (z < _kernelRadius) || (z > this->_input->GetZ() - 1 - _kernelRadius)) {
            this->_output->Put(x, y, z, t, this->_input->Get(x, y, z, t));
          } else {
        	  vector<VoxelType> voxels;
        	  for( int zz = z - _kernelRadius; zz <= z + _kernelRadius; ++zz )
        	  {
        		  for( int yy = y - _kernelRadius; yy <= y + _kernelRadius; ++yy )
        		  {
        			  for( int xx = x - _kernelRadius; xx <= x + _kernelRadius; ++xx )
        			  {
        				  if( this->_mask->Get(xx,yy,zz,t) ) {
        					  voxels.push_back( this->_input->Get(xx,yy,zz,t) );
        				  }
        			  }
        		  }
        	  }
        	  if( voxels.size() )
        	  {
        		  sort(voxels.begin(), voxels.end());
        		  this->_output->Put(x, y, z, t, voxels[voxels.size()/2]);
        	  }
          }
        }
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();
}

template class irtkMedianFilter<irtkBytePixel>;
template class irtkMedianFilter<irtkGreyPixel>;
template class irtkMedianFilter<irtkRealPixel>;
