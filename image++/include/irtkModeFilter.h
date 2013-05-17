
#ifndef _IRTKMODEFILTER_H

#define _IRTKMODEFILTER_H

#include <irtkImageToImage.h>

/**
 * Class for applying mode filter to what should be label images.
 *
 * Assign to each voxel, the modal label of those within a neighbourhood.
 */

template <class VoxelType> class irtkModeFilter : public irtkImageToImage<VoxelType>
{

protected:

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Initialize the filter
  virtual void Initialize();

  /// What connectivity to assume when running the filter.
  irtkConnectivityType _Connectivity;

  // List of voxel offsets of the neighbourhood.
  irtkNeighbourhoodOffsets _offsets;

public:

  /// Constructor
  irtkModeFilter();

  /// Destructor
  ~irtkModeFilter();

  /// Run mode filter
  virtual void Run();

  SetMacro(Connectivity, irtkConnectivityType);

  GetMacro(Connectivity, irtkConnectivityType);
};

#endif
