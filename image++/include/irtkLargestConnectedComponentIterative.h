/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKLARGESTCONNECTEDCOMPONENTITERATIVE_H

#define _IRTKLARGESTCONNECTEDCOMPONENTITERATIVE_H

#include <irtkImageToImage.h>

/**
 * Class for extracting the largest connected component from a labelled image
 *
 * This class defines and implements the extraction of the largest
 * connected component from a labelled image.  Uses an iterative method so
 * takes longer than the class irtkLargestConnectedComponent.  The stack
 * size needed, however, is lower.
 *
 */

template <class VoxelType> class irtkLargestConnectedComponentIterative : public irtkImageToImage<VoxelType>
{

  /// Label used to identify labels of interest
  VoxelType _TargetLabel;

  /// Size of largest cluster with the target label.
  int _largestClusterSize;

  // What label is used for the largest cluster during computation.
  VoxelType _largestClusterLabel;

  // Record of all cluster sizes.
  int *_ClusterSizes;

  // Number of clusters that match the target label.
  int _NumberOfClusters;

  /// Modes
  bool _Mode2D;

  bool _AllClustersMode;

protected:

  /// Run filter
  virtual void Run3D();
  virtual void Run2D();

  // Helper functions
  virtual void ResetMarks();

  virtual void SelectLargestCluster();

  int CheckAdjacency2D(VoxelType& markA, VoxelType& markB);
  int CheckAdjacency3D(VoxelType& markA, VoxelType& markB);

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Returns the name of the class
  virtual const char *NameOfClass();

public:

  /// Constructor
  irtkLargestConnectedComponentIterative(VoxelType = 0);

  /// Destructor
  ~irtkLargestConnectedComponentIterative();

  /// Set label sought
  SetMacro(TargetLabel, VoxelType);

  /// Get label sought
  GetMacro(TargetLabel, VoxelType);

  /// Set mode
  SetMacro(Mode2D, bool);

  /// Get mode
  GetMacro(Mode2D, bool);

  /// Set mode
  SetMacro(AllClustersMode, bool);

  /// Get mode
  GetMacro(AllClustersMode, bool);

  /// Run filter
  virtual void Run();
};

#endif
