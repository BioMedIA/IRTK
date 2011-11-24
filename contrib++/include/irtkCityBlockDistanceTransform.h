/*=========================================================================

=========================================================================*/

/*=========================================================================

 Find the City Block (Manhattan, L1) distance for all object voxels in an
 image from the boundary. Object voxels have a value greater than zero and
 background voxels are the rest.  The distance map is initialised to zero.
 In each iteration, the distance map is incremented by 1 for all object
 voxels. The border voxels are removed and the the next iteration starts.
 If any dimension is a singleton, a 2D version is applied.
 PA, 2011 07 04

=========================================================================*/

#ifndef _IRTKCITYBLOCKDISTANCETRANSFORM_H

#define _IRTKCITYBLOCKDISTANCETRANSFORM_H

#include <irtkImageToImage.h>

/// In the 2D case, a flip may be necessary so that the singleton dimension
/// is the z-direction. This makes processing easier.
typedef enum { irtkFlipNone, irtkFlipXY,  irtkFlipXZ,  irtkFlipYZ} irtkFlipType;

template <class VoxelType> class irtkCityBlockDistanceTransform : public irtkImageToImage<VoxelType>
{

protected:

	/// Storage for the object voxels that can be updated.
	irtkGreyImage *_data;

  /// List of voxel offsets of the 6 neighbourhood of a voxel. Enables checking to
	/// see if there is a face-neighbour labelled background during each iteration.
  irtkNeighbourhoodOffsets _offsets;

  irtkFlipType _flipType;

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Requires buffering
  virtual bool RequiresBuffering();

  /// Run distance transform in 2D
  virtual void Run2D();

  /// Run distance transform in 3D
  virtual void Run3D();

  /// Initialize the filter
  virtual void Initialize();

  /// Initialize the filter
  virtual void Initialize2D();

  /// Initialize the filter
  virtual void Initialize3D();

  /// Finalize the filter
  virtual void Finalize();

public:

  /// Default constructor
  irtkCityBlockDistanceTransform();

  /// Destructor (empty).
  ~irtkCityBlockDistanceTransform() {};

  /// Run distance transform
  virtual void Run();

};

template <class VoxelType> inline const char *irtkCityBlockDistanceTransform<VoxelType>::NameOfClass()
{
  return "irtkCityBlockDistanceTransform";
}

template <class VoxelType> inline bool irtkCityBlockDistanceTransform<VoxelType>::RequiresBuffering()
{
  return false;
}

#endif
