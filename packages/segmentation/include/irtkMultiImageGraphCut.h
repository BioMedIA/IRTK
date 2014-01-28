/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#ifndef _irtkMultiImageGraphCut_H

#define _irtkMultiImageGraphCut_H

#include <irtkImage.h>
#include <GCoptimization.h>
/**
 * Abstract base class for any general image to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image as output. Each
 * derived class has to implement all abstract member functions.
 */

template <class VoxelType> class irtkMultiImageGraphCut : public irtkObject
{

protected:

  /// Number of input images
  int _numberOfImages;

  /// number of labels
  int _labels;

  double *_datacost;

  /// spacing
  double *_dx,*_dy,*_dz,_dt;

  /// Number of Voxels
  int _totalVoxel;

  /// Image node offset
  int *_imageoffset;

  /// Input image for filter
  irtkGenericImage<VoxelType> **_input;

  /// Output image for filter
  irtkGenericImage<VoxelType> **_output;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /** Finalize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Finalize();

  /// add regionweight to graph
  virtual void AddBoundaryTerm(GCoptimizationGeneralGraph *graph, int count, 
	  int i,int j, int k, int l, int n,
	  int xoff, int yoff, int zoff, int toff, double divide);

  /// add multiimage weight to graph
  virtual void AddImageTerm(GCoptimizationGeneralGraph *graph, int count, 
	  int count2, double divide);

  /// geometry mode 0 no connection 1 x 2 xy 3 xyz 4 xyzt
  int _mode;

public:

  /// Constructor
  irtkMultiImageGraphCut();

  /// Deconstuctor
  virtual ~irtkMultiImageGraphCut();

  /// Set input image for filter
  virtual void SetInput (int, irtkGenericImage<VoxelType> **,int, double *);

  /// Set output image for filter
  virtual void SetOutput(irtkGenericImage<VoxelType> **);

  /// Set Geometry Mode
  virtual void SetMode(int);

  /// Set dt
  SetMacro(dt,double);

  /// Run filter on entire image lambda is the coefficent for region term.
  virtual void Run(double = 1, int = 0);

  /// Get Min Value
  virtual VoxelType MinValue(VoxelType,VoxelType);

    /// Get Max Value
  virtual VoxelType MaxValue(VoxelType,VoxelType);

  /// Returns the name of the class
  virtual const char *NameOfClass();
};

template <class VoxelType> const char *irtkMultiImageGraphCut<VoxelType>::NameOfClass()
{
  return "irtkMultiImageGraphCut";
}

template <class VoxelType> void irtkMultiImageGraphCut<VoxelType>::SetMode(int i)
{
  _mode = i;
}

template <class VoxelType> VoxelType irtkMultiImageGraphCut<VoxelType>::MinValue(VoxelType first, VoxelType second)
{
  if(first > second)
	  return second;
  else
	  return first;
}

template <class VoxelType> VoxelType irtkMultiImageGraphCut<VoxelType>::MaxValue(VoxelType first, VoxelType second)
{
  if(first < second)
	  return second;
  else
	  return first;
}

#endif
