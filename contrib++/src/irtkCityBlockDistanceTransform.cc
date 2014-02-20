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

#include <irtkCityBlockDistanceTransform.h>

template <class VoxelType> irtkCityBlockDistanceTransform<VoxelType>::irtkCityBlockDistanceTransform() : irtkImageToImage<VoxelType>()
{
	this->_flipType = irtkFlipNone;
}

template <class VoxelType> void irtkCityBlockDistanceTransform<VoxelType>::Initialize2D()
{
	int i, j, l, nx, ny, nz, nt;
	irtkImageAttributes attr;

	// Make an array with an extra layer in each direction.
	attr = this->_input->GetImageAttributes();

	// One of the dimensions should be a singleton.
	if (attr._x != 1 && attr._y != 1 && attr._z != 1){
		cerr << "irtkCityBlockDistanceTransform<VoxelType>::Initialize2D : One spatial dimension should be a singleton." << endl;
		exit(1);
	}

	// Ensure the singleton spatial dimension is the z-dimension.

	// Default.
	this->_flipType = irtkFlipNone;

	if (attr._x == 1){
		this->_input->FlipXZ(0);
		this->_output->FlipXZ(0);
		this->_flipType = irtkFlipXZ;
	} else if (attr._y == 1) {
		this->_input->FlipYZ(0);
		this->_output->FlipYZ(0);
		this->_flipType = irtkFlipYZ;
	}


	// Get the attributes again as they may have changed.
	attr = this->_input->GetImageAttributes();

	// Increment the planar dimensions.
	attr._x += 2;
	attr._y += 2;

  this->_data = new irtkGreyImage(attr);

  nx = attr._x;
  ny = attr._y;
	nz = attr._z; // == 1
	nt = attr._t;

  for (l = 0; l < nt; ++l){

  	// Outer rows and columns are duplicates of outer rows and
  	// columns of input image.
  	for (j = 1; j < ny - 1; ++j){
  		this->_data->Put(0   , j, 0, l,
  				this->_input->Get(0   , j-1, 0, l) > 0 ? 1 : 0);
  		this->_data->Put(nx-1, j, 0, l,
  				this->_input->Get(nx-3, j-1, 0, l) > 0 ? 1 : 0);
  	}

  	for (i = 1; i < nx - 1; ++i){
  		this->_data->Put(i, 0   , 0, l,
  				this->_input->Get(i-1, 0   , 0, l) > 0 ? 1 : 0);
  		this->_data->Put(i, ny-1, 0, l,
  				this->_input->Get(i-1, ny-3, 0, l) > 0 ? 1 : 0);
  	}

  	// Copy original image into interior.
  	for (j = 1; j < ny - 1; ++j){
  		for (i = 1; i < nx - 1; ++i){
  			this->_data->Put(i, j, 0, l,
  					this->_input->Get(i-1, j-1, 0, l) > 0 ? 1 : 0);
  		}
  	}
  }

  // Only need to use the first four offsets in the 2D case.
	this->_offsets.Initialize(nx, ny, CONNECTIVITY_06);
}


template <class VoxelType> void irtkCityBlockDistanceTransform<VoxelType>::Initialize3D()
{
	int i, j, k, l, nx, ny, nz, nt;
	irtkImageAttributes attr;

	// Make an array with an extra layer in each direction.
	attr = this->_input->GetImageAttributes();
	attr._x += 2;
	attr._y += 2;
	attr._z += 2;

  this->_data = new irtkGreyImage(attr);

  nx = this->_data->GetX();
	ny = this->_data->GetY();
	nz = this->_data->GetZ();
	nt = this->_data->GetT();

  for (l = 0; l < nt; ++l){

  	// Outer layers are duplicates of outer layer of input image.
  	for (j = 1; j < ny - 1; ++j){
  		for (i = 1; i < nx - 1; ++i){
  			this->_data->Put(i, j, 0   , l,
  					this->_input->Get(i-1, j-1, 0   , l) > 0 ? 1 : 0);
  			this->_data->Put(i, j, nz-1, l,
  					this->_input->Get(i-1, j-1, nz-3, l) > 0 ? 1 : 0);
  		}
  	}

  	for (k = 1; k < nz - 1; ++k){
  		for (i = 1; i < nx - 1; ++i){
  			this->_data->Put(i, 0   , k, l,
  					this->_input->Get(i-1, 0   , k-1, l) > 0 ? 1 : 0);
  			this->_data->Put(i, ny-1, k, l,
  					this->_input->Get(i-1, ny-3, k-1, l) > 0 ? 1 : 0);
  		}
  	}

  	for (k = 1; k < nz - 1; ++k){
  		for (j = 1; j < ny - 1; ++j){
  			this->_data->Put(0   , j, k, l,
  					this->_input->Get(0   , j-1, k-1, l) > 0 ? 1 : 0);
  			this->_data->Put(nx-1, j, k, l,
  					this->_input->Get(nx-3, j-1, k-1, l) > 0 ? 1 : 0);
  		}
  	}

  	// Copy original image into interior.
  	for (k = 1; k < nz - 1; ++k){
  		for (j = 1; j < ny - 1; ++j){
    		for (i = 1; i < nx - 1; ++i){
    			this->_data->Put(i, j, k, l,
    					this->_input->Get(i-1, j-1, k-1, l) > 0 ? 1 : 0);
    		}
  		}
  	}


  }

  // Offsets for searching in the neighbourhood of a voxel.
	this->_offsets.Initialize(nx, ny, CONNECTIVITY_06);
}

template <class VoxelType> void irtkCityBlockDistanceTransform<VoxelType>::Initialize()
{
	int nx, ny, nz;

  // Do the initial set up
  this->irtkImageToImage<VoxelType>::Initialize();

	nx = this->_input->GetX();
	ny = this->_input->GetY();
	nz = this->_input->GetZ();

	if (nx == 1 || ny == 1 || nz == 1){
		this->Initialize2D();
	} else {
		this->Initialize3D();
	}

}

template <class VoxelType> void irtkCityBlockDistanceTransform<VoxelType>::Finalize()
{
	// If image is 2D, may need to undo a dimension flip.
	if (this->_data->GetZ() == 1){
		if (this->_flipType == irtkFlipXZ){
			this->_output->FlipXZ(0);
		} else if (this->_flipType == irtkFlipYZ) {
			this->_output->FlipYZ(0);
		}
	}

	delete this->_data;

	this->irtkImageToImage<VoxelType>::Finalize();
}

template <class VoxelType> void irtkCityBlockDistanceTransform<VoxelType>::Run2D()
{
	int nx, ny, nz, nt;
	int i, j, l, m, n;
	int borderVoxelCount = 0, objectVoxelCount = 0;
	int *borderIndices;
	VoxelType val;

  irtkGreyPixel *ptr2current, *ptr2offset, *ptr;

	nx = this->_data->GetX();
	ny = this->_data->GetY();
	nz = this->_data->GetZ();
	nt = this->_data->GetT();

	if (nz != 1){
		cerr << "irtkCityBlockDistanceTransform<VoxelType>::Run2D() : Expect one slice only" << endl;
		exit(1);
	}

	// Storage for recording the indices of border voxels.
	borderIndices = new int[nx * ny];

	for (l = 0; l < nt; ++l){
		do{ // while (objectVoxelCount > 0);

			// Increment the distance for all current object voxels.
			objectVoxelCount = 0;

			for (j = 1; j < ny - 1; ++j){
				for (i = 1; i < nx - 1; ++i){
					if (this->_data->Get(i, j, 0, l) > 0){
						++objectVoxelCount;
						val = this->_output->Get(i-1, j-1, 0, l);
						this->_output->Put(i-1, j-1, 0, l, 1 + val);
					}
				}
			}

			if (objectVoxelCount == 0){
				cerr << "irtkCityBlockDistanceTransform<VoxelType>::Run2D() : No object voxels." << endl;
				break;
			}

			// Remove the border from the current object(s).
			borderVoxelCount = 0;
			ptr = this->_data->GetPointerToVoxels(0, 0, 0, l);

			for (j = 1; j < ny - 1; ++j){
				for (i = 1; i < nx - 1; ++i){

					if (this->_data->Get(i, j, 0, l) > 0){

						// Object voxel. Check face neighbourhood for background.
						ptr2current = this->_data->GetPointerToVoxels(i, j, 0, l);

						// 2D case: Only need to check first 4 neighbours in offset list.
						for (m = 0; m < 4; ++m){

							ptr2offset = ptr2current + this->_offsets(m);
							if (*ptr2offset < 1){
								// Border voxel.
								borderIndices[borderVoxelCount] = ptr2current - ptr;

								++borderVoxelCount;
								break;
							}
						}

					}
				}
			}

			// Remove the current border voxels.
			ptr = this->_data->GetPointerToVoxels(0, 0, 0, l);

			for (n = 0; n < borderVoxelCount; ++n){
				ptr2current = ptr + borderIndices[n];
				(*ptr2current) = 0;
			}

			// Update count of object voxels.
			objectVoxelCount -= borderVoxelCount;

		} while (objectVoxelCount > 0);

	}

	delete [] borderIndices;
}

template <class VoxelType> void irtkCityBlockDistanceTransform<VoxelType>::Run3D()
{
	int nx, ny, nz, nt;
	int i, j, k, l, m, n;
	int borderVoxelCount = 0, objectVoxelCount = 0;
	int *borderIndices;
	VoxelType val;

  irtkGreyPixel *ptr2current, *ptr2offset, *ptr2start;

	nx = this->_data->GetX();
	ny = this->_data->GetY();
	nz = this->_data->GetZ();
	nt = this->_data->GetT();

	// Storage for recording the indices of border voxels.
	borderIndices = new int[nx * ny * nz];

	for (l = 0; l < nt; ++l){


		do{ // while (objectVoxelCount > 0);

			// Increment the distance for all current object voxels.
			objectVoxelCount = 0;
			for (k = 1; k < nz - 1; ++k){
				for (j = 1; j < ny - 1; ++j){
					for (i = 1; i < nx - 1; ++i){
						if (this->_data->Get(i, j, k, l) > 0){
							++objectVoxelCount;
							val = this->_output->Get(i-1, j-1, k-1, l);
							this->_output->Put(i-1, j-1, k-1, l, 1 + val);
						}

					}
				}
			}

			if (objectVoxelCount == 0){
				cerr << "irtkCityBlockDistanceTransform<VoxelType>::Run3D() : No object voxels." << endl;
				break;
			}

			// Remove the border from the current object.
			borderVoxelCount = 0;
			ptr2start = this->_data->GetPointerToVoxels(0, 0, 0, l);

			for (k = 1; k < nz - 1; ++k){
				for (j = 1; j < ny - 1; ++j){
					for (i = 1; i < nx - 1; ++i){
						if (this->_data->Get(i, j, k, l) > 0){

							// Object voxel. Check face neighbourhood for background (6-neighbourhood).
							ptr2current = this->_data->GetPointerToVoxels(i, j, k, l);

							for (m = 0; m < 6; ++m){
								ptr2offset = ptr2current + this->_offsets(m);
								if (*ptr2offset < 1){
									// Border voxel.
									borderIndices[borderVoxelCount] = ptr2current - ptr2start;

									++borderVoxelCount;
									break;
								}
							}
						}

					}
				}
			}


			// Remove the border voxels.
			for (n = 0; n < borderVoxelCount; ++n){
				ptr2current = ptr2start + borderIndices[n];
				(*ptr2current) = 0;
			}

			// Update count of object voxels.
			objectVoxelCount -= borderVoxelCount;


		} while (objectVoxelCount > 0);

	}

	delete [] borderIndices;
}

template <class VoxelType> void irtkCityBlockDistanceTransform<VoxelType>::Run()
{
	int nx, ny, nz;

	// Do the initial set up
	this->Initialize();

	// Calculate image dimensions
	nx = this->_input->GetX();
	ny = this->_input->GetY();
	nz = this->_input->GetZ();

	if (nx == 1 || ny == 1 || nz == 1) {
		// Calculate 2D distance transform
		this->Run2D();
	} else {
		this->Run3D();
	}


  // Do the final cleaning up
  this->Finalize();
}


template class irtkCityBlockDistanceTransform<irtkRealPixel>;

template class irtkCityBlockDistanceTransform<irtkGreyPixel>;
