/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKLARGESTCONNECTEDCOMPONENT_H

#define _IRTKLARGESTCONNECTEDCOMPONENT_H

#include <irtkImageToImage.h>

/**
 * Class for extracting the largest connected component from a labelled image
 *
 * This class defines and implements the extraction of the largest connected component
 * from a labelled image.
 *
 */

template <class VoxelType> class irtkLargestConnectedComponent : public irtkImageToImage<VoxelType>
{

  /// Size of current cluster
  int _currentClusterSize;

  /// Size of largest cluster
  int _largestClusterSize;

  /// Label used to identify labels of interest
  VoxelType _ClusterLabel;

  /// Mode
  bool _Mode2D;

protected:

  /// Recursive region growing
  void Grow2D(int, int, int, int);

  /// Recursive region growing
  void Grow3D(int, int, int, int);

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor
  irtkLargestConnectedComponent(VoxelType = 0);

  /// Destructor
  ~irtkLargestConnectedComponent();

  /// Set sigma
  SetMacro(ClusterLabel, VoxelType);

  /// Get sigma
  GetMacro(ClusterLabel, VoxelType);

  /// Set mode
  SetMacro(Mode2D, bool);

  /// Get mode
  GetMacro(Mode2D, bool);

  /// Run filter
  virtual void Run();

};

#endif
