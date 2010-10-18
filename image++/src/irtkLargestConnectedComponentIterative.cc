/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <map>

#include <irtkImage.h>

#include <irtkLargestConnectedComponentIterative.h>

typedef map<int, int> countMap;

// Constructor.
template <class VoxelType> irtkLargestConnectedComponentIterative<VoxelType>::irtkLargestConnectedComponentIterative(VoxelType TargetLabel)
{
  _largestClusterSize  = 0;
  _largestClusterLabel = 0;
  _AllClustersMode     = false;
  _Mode2D              = false;
  _TargetLabel         = TargetLabel;
  _NumberOfClusters    = 0;
  _ClusterSizes        = NULL;
}

template <class VoxelType> irtkLargestConnectedComponentIterative<VoxelType>::~irtkLargestConnectedComponentIterative(void)
{
  delete [] _ClusterSizes;
}

template <class VoxelType> bool irtkLargestConnectedComponentIterative<VoxelType>::RequiresBuffering(void)
{
  return true;
}

template <class VoxelType> const char *irtkLargestConnectedComponentIterative<VoxelType>::NameOfClass()
{
  return "irtkLargestConnectedComponentIterative";
}

template <class VoxelType> int irtkLargestConnectedComponentIterative<VoxelType>::CheckAdjacency2D(VoxelType& markA, VoxelType& markB)
{
  int foundAdjacency = false;
  int i, j, k = 0, l = 0;
  int xdim, ydim;

  xdim = this->_input->GetX();
  ydim = this->_input->GetY();

  VoxelType temp;

  for (j = 1; j < ydim - 1; ++j) {
    for (i = 1; i < xdim - 1; ++i) {
      markA = this->_output->Get(i, j, k, l);

      if (markA < 1)
        continue;

      // Check the 4 neighbourhood.
      temp = this->_output->Get(i + 1, j, k, l);
      if (temp > 0 && markA != temp && foundAdjacency == false) {
        foundAdjacency = true;
        markB = temp;
      }

      temp = this->_output->Get(i - 1, j, k, l);
      if (temp > 0 && markA != temp && foundAdjacency == false) {
        foundAdjacency = true;
        markB = temp;
      }

      temp = this->_output->Get(i, j + 1, k, l);
      if (temp > 0 && markA != temp && foundAdjacency == false) {
        foundAdjacency = true;
        markB = temp;
      }

      temp = this->_output->Get(i, j - 1, k, l);
      if (temp > 0 && markA != temp && foundAdjacency == false) {
        foundAdjacency = true;
        markB = temp;
      }

      if (foundAdjacency == true)
        break;
    }
    if (foundAdjacency == true)
      break;
  }

  return foundAdjacency;
}

template <class VoxelType> int irtkLargestConnectedComponentIterative<VoxelType>::CheckAdjacency3D(VoxelType& markA, VoxelType& markB)
{
  int foundAdjacency = false;
  int i, j, k, l = 0;

  int xdim, ydim, zdim;
  xdim = this->_input->GetX();
  ydim = this->_input->GetY();
  zdim = this->_input->GetZ();
  VoxelType temp;

  for (k = 1; k < zdim - 1; ++k) {
    for (j = 1; j < ydim - 1; ++j) {
      for (i = 1; i < xdim - 1; ++i) {
        markA = this->_output->Get(i, j, k, l);

        if (markA < 1)
          continue;

        // Check the 6 neighbourhood.
        temp = this->_output->Get(i + 1, j, k, l);
        if (temp > 0 && markA != temp && foundAdjacency == false) {
          foundAdjacency = true;
          markB = temp;
        }

        temp = this->_output->Get(i - 1, j, k, l);
        if (temp > 0 && markA != temp && foundAdjacency == false) {
          foundAdjacency = true;
          markB = temp;
        }

        temp = this->_output->Get(i, j + 1, k, l);
        if (temp > 0 && markA != temp && foundAdjacency == false) {
          foundAdjacency = true;
          markB = temp;
        }

        temp = this->_output->Get(i, j - 1, k, l);
        if (temp > 0 && markA != temp && foundAdjacency == false) {
          foundAdjacency = true;
          markB = temp;
        }

        temp = this->_output->Get(i, j, k + 1, l);
        if (temp > 0 && markA != temp && foundAdjacency == false) {
          foundAdjacency = true;
          markB = temp;
        }

        temp = this->_output->Get(i, j, k - 1, l);
        if (temp > 0 && markA != temp && foundAdjacency == false) {
          foundAdjacency = true;
          markB = temp;
        }

        if (foundAdjacency == true)
          break;
      }
      if (foundAdjacency == true)
        break;
    }
    if (foundAdjacency == true)
      break;
  }

  return foundAdjacency;
}

template <class VoxelType> void irtkLargestConnectedComponentIterative<VoxelType>::SelectLargestCluster()
{
  int i, voxels;
  VoxelType *ptr;
  voxels = this->_output->GetNumberOfVoxels();
  ptr = this->_output->GetPointerToVoxels();

  for (i = 0; i < voxels; ++i) {

    if (*ptr == _largestClusterLabel) {
      *ptr = 1;
    } else {
      *ptr = 0;
    }

    ++ptr;
  }
}

// Replace the voxels labeled with a variety of marks with a sequence starting from 1.
template <class VoxelType> void irtkLargestConnectedComponentIterative<VoxelType>::ResetMarks()
{
  int voxels, i, j;
  VoxelType *ptr;
  int *oldLabels;

  countMap counts;
  map<int, int>::iterator iter;

  voxels = this->_output->GetNumberOfVoxels();
  ptr = this->_output->GetPointerToVoxels();

  // Count up the number of labels for each positive label.
  for (i = 0; i < voxels; ++i) {
    if (*ptr > 0) {
      ++counts[*ptr];
    }
    ++ptr;
  }

  _NumberOfClusters = counts.size();

  if (_NumberOfClusters < 1) {
    cerr << "ResetMarks : There are no clusters." << endl;
    exit(1);
  }

  oldLabels = new int[_NumberOfClusters];
  _ClusterSizes = new int[_NumberOfClusters];

  cout << "There are " << _NumberOfClusters << " clusters." << endl;

  // What labels are used and which has the largest cluster?
  _largestClusterSize = 0;
  i = 0;

  for (iter = counts.begin(); iter != counts.end(); ++iter) {

    oldLabels[i] = iter->first;
    _ClusterSizes[i] = iter->second;
    ++i;

    if (iter->second > _largestClusterSize) {
      _largestClusterSize  = iter->second;
      _largestClusterLabel = i;
    }
  }

  // Now replace the labels with a sequence 1, 2, ...
  ptr = this->_output->GetPointerToVoxels();

  for (i = 0; i < voxels; ++i) {
    if (*ptr > 0) {
      for (j = 0; j < _NumberOfClusters; ++j) {
        if (*ptr == oldLabels[j]) {
          *ptr = j + 1;
          break;
        }
      }
    }
    ++ptr;
  }

  delete [] oldLabels;
}

template <class VoxelType> void irtkLargestConnectedComponentIterative<VoxelType>::Run2D()
{
  int i, j, k = 0, l = 0;
  int xdim, ydim;
  int voxelsMarked = true;
  VoxelType currentMark = 0;
  VoxelType markA, markB;
  int foundAdjacency = true;
  VoxelType in, out;

  xdim = this->_input->GetX();
  ydim = this->_input->GetY();

  // Iterate over the input image, placing marks in the matching output
  // images where the input has the required label and an adjacent,
  // previously visited, voxel has the same mark.
  int marksAssigned;

  while (voxelsMarked == true) {
    voxelsMarked = false;
    currentMark++;

    marksAssigned = 0;

    for (j = 0; j < ydim; ++j) {
      for (i = 0; i < xdim; ++i) {

        in  = this->_input->Get(i, j, k, l);
        out = this->_output->Get(i, j, k, l);

        // Have we found an unmarked target label?
        if (in == _TargetLabel && out == 0) {

          if (voxelsMarked == false) {
            // This is the first one found.
            this->_output->Put(i, j, k, l, currentMark);
            voxelsMarked = true;
            ++marksAssigned;
          } else if ( (i > 0 && this->_output->Get(i - 1, j, k, l) == currentMark) ||
                      (j > 0 && this->_output->Get(i, j - 1, k, l) == currentMark) ) {
            this->_output->Put(i, j, k, l, currentMark);
            ++marksAssigned;
          }
        }
      }
    }
  }

  // Now relabel adjacent marks to a common mark.
  while (foundAdjacency == true) {
    foundAdjacency = CheckAdjacency2D(markA, markB);

    if (foundAdjacency == true) {
      // Update both marks to the next available mark.
      currentMark++;

      for (j = 0; j < ydim; ++j) {
        for (i = 0; i < xdim; ++i) {
          if (this->_output->Get(i, j, k, l) == markA ||
              this->_output->Get(i, j, k, l) == markB) {
            this->_output->Put(i, j, k, l, currentMark);
          }
        }
      }
    }
  }

  this->ResetMarks();

  if (this->_AllClustersMode == false) {
    // We want only the largest cluster.
    this->SelectLargestCluster();
  }
}


template <class VoxelType> void irtkLargestConnectedComponentIterative<VoxelType>::Run3D()
{
  int i, j, k, l = 0;
  int xdim, ydim, zdim;
  int voxelsMarked = true;
  VoxelType currentMark = 0;
  VoxelType markA, markB;
  int foundAdjacency = true;
  VoxelType in, out;

  xdim = this->_input->GetX();
  ydim = this->_input->GetY();
  zdim = this->_input->GetZ();

  // Iterate over the input image, placing marks in the matching output
  // images where the input has the required label and an adjacent,
  // previously visited, voxel has the same mark.
  int marksAssigned;

  while (voxelsMarked == true) {
    voxelsMarked = false;
    currentMark++;

    marksAssigned = 0;

    for (k = 0; k < zdim; ++k) {
      for (j = 0; j < ydim; ++j) {
        for (i = 0; i < xdim; ++i) {

          in  = this->_input->Get(i, j, k, l);
          out = this->_output->Get(i, j, k, l);

          // Have we found an unmarked target label?
          if (in == _TargetLabel && out == 0) {

            if (voxelsMarked == false) {
              // This is the first one found.
              this->_output->Put(i, j, k, l, currentMark);
              voxelsMarked = true;
              ++marksAssigned;
            } else if ( (i > 0 && this->_output->Get(i - 1, j, k, l) == currentMark) ||
                        (j > 0 && this->_output->Get(i, j - 1, k, l) == currentMark) ||
                        (k > 0 && this->_output->Get(i, j, k - 1, l) == currentMark) ) {
              this->_output->Put(i, j, k, l, currentMark);
              ++marksAssigned;
            }
          }
        }
      }
    }
  }

  // Now relabel adjacent marks to a common mark.
  while (foundAdjacency == true) {
    foundAdjacency = CheckAdjacency3D(markA, markB);

    if (foundAdjacency == true) {
      // Update both marks to the next available mark.
      currentMark++;

      for (k = 0; k < zdim; ++k) {
        for (j = 0; j < ydim; ++j) {
          for (i = 0; i < xdim; ++i) {
            if (this->_output->Get(i, j, k, l) == markA ||
                this->_output->Get(i, j, k, l) == markB) {
              this->_output->Put(i, j, k, l, currentMark);
            }
          }
        }
      }
    }
  }

  this->ResetMarks();

  if (this->_AllClustersMode == false) {
    // We want only the largest cluster.
    this->SelectLargestCluster();
  }
}

template <class VoxelType> void irtkLargestConnectedComponentIterative<VoxelType>::Run()
{
  int x, y, z, t = 0;

  // Do the initial set up
  this->Initialize();

  // Check if the input is 4D
  if (this->_input->GetT() > 1) {
    cerr << "irtkLargestConnectedComponentIterative<VoxelType>::Run(): 4D images not yet supported\n" << endl;
    exit(1);
  }

  // Reset the output.
  for (z = 0; z < this->_input->GetZ(); z++) {
    for (y = 0; y < this->_input->GetY(); y++) {
      for (x = 0; x < this->_input->GetX(); x++) {
        this->_output->Put(x, y, z, t, 0);
      }
    }
  }

  // Main calls.
  if (this->_Mode2D == true) {
    if (this->_input->GetZ() != 1) {
      cerr << "irtkLargestConnectedComponentIterative::Run : ";
      cerr << "2D mode selected but image has more than one slice in the z direction." << endl;
      exit(1);
    }
    this->Run2D();
  } else {
    this->Run3D();
  }

  // Do the final cleaning up
  this->Finalize();
}

template class irtkLargestConnectedComponentIterative<irtkBytePixel>;
template class irtkLargestConnectedComponentIterative<irtkGreyPixel>;
//template class irtkLargestConnectedComponentIterative<irtkRealPixel>;
