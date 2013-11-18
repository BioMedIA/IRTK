/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date: 2013-01-14$
Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>
#include <irtkMeanFilter.h>

template <class VoxelType> irtkMeanFilter<VoxelType>::irtkMeanFilter()
{
	// Default kernel radius.
	this->_radius_x = 2;
	this->_radius_y = 2;
	this->_radius_z = 2;
}

template <class VoxelType> irtkMeanFilter<VoxelType>::~irtkMeanFilter(void)
{
}

template <class VoxelType> const char *irtkMeanFilter<VoxelType>::NameOfClass()
{
	return "irtkMeanFilter";
}

template <class VoxelType> void irtkMeanFilter<VoxelType>::Initialize()
{
	// Do the initial set up
	this->irtkImageToImage<VoxelType>::Initialize();
}

template <class VoxelType> bool irtkMeanFilter<VoxelType>::RequiresBuffering(void)
{
	return true;
}

template <class VoxelType> void irtkMeanFilter<VoxelType>::Run()
{
	int x, y, z, t;

	//int dim = 2*_kernelRadius + 1;

	// Do the initial set up
	this->Initialize();

	for (t = 0; t < this->_input->GetT(); t++) {
		for (z = 0; z < this->_input->GetZ(); z++) {
			for (y = 0; y < this->_input->GetY(); y++) {
				for (x = 0; x < this->_input->GetX(); x++) {
					double sum = 0, count = 0;
					for( int zz = z - _radius_z; zz <= z + _radius_z; ++zz )
					{
						for( int yy = y - _radius_y; yy <= y + _radius_y; ++yy )
						{
							for( int xx = x - _radius_x; xx <= x + _radius_x; ++xx )
							{
								if ((xx < 0) || (xx > this->_input->GetX() - 1) ||
									(yy < 0) || (yy > this->_input->GetY() - 1) ||
									(zz < 0) || (zz > this->_input->GetZ() - 1)) {
								} else {
									sum += this->_input->Get(xx, yy, zz, t);
									count += 1.0;
								}
							}
						}
					}
					if( count > 0 )
					{
						this->_output->Put(x, y, z, t, sum/count);
					}else
					{
						this->_output->Put(x, y, z, t, this->_input->Get(x, y, z, t));				  
					}
				}
			}
		}
	}

	// Do the final cleaning up
	this->Finalize();
}

template <class VoxelType> void irtkMeanFilter<VoxelType>::SetkernelRadius(double radius)
{
	if(_input != NULL){
		this->_radius_x = round(radius / _input->GetXSize());
		this->_radius_y = round(radius / _input->GetYSize());
		this->_radius_z = round(radius / _input->GetZSize());
	}else{
		this->_radius_x = radius;
		this->_radius_y = radius;
		this->_radius_z = radius;
	}

	if(this->_radius_x < 1) this->_radius_x = 1;
	if(this->_radius_y < 1) this->_radius_y = 1;
	if(this->_radius_z < 1) this->_radius_z = 1;
}

template class irtkMeanFilter<irtkBytePixel>;
template class irtkMeanFilter<irtkGreyPixel>;
template class irtkMeanFilter<irtkRealPixel>;
