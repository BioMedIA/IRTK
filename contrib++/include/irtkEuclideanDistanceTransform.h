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

#ifndef _IRTKEUCLIDEANDISTANCETRANSFORM_H

#define _IRTKEUCLIDEANDISTANCETRANSFORM_H

#define EDT_MAX_IMAGE_DIMENSION 26754
#define EDT_MAX_DISTANCE_SQUARED 2147329548
#define EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC 2147329548

#include <irtkImageToImage.h>

template <class VoxelType> class irtkEuclideanDistanceTransform : public irtkImageToImage<VoxelType>
{

public:

  /// 2D or 3D distance transform
  enum irtkDistanceTransformMode { irtkDistanceTransform2D, irtkDistanceTransform3D };

protected:

  /// 2D or 3D distance transform
  irtkDistanceTransformMode _distanceTransformMode;

  /// Calculate the Vornoi diagram
  int edtVornoiEDT(long *, long);

  /// Calculate 2D distance transform
  void edtComputeEDT_2D(char *, long *, long, long);

  /// Calculate 3D distance transform
  void edtComputeEDT_3D(char *, long *, long, long, long);

  /// Calculate the Vornoi diagram for anisotripic voxel sizes
  int edtVornoiEDT_anisotropic(irtkRealPixel *, long, double);

  /// Calculate 2D distance transform for anisotripic voxel sizes
  void edtComputeEDT_2D_anisotropic(irtkRealPixel *, irtkRealPixel *, long, long, double, double);

  /// Calculate 3D distance transform for anisotripic voxel sizes
  void edtComputeEDT_3D_anisotropic(irtkRealPixel *, irtkRealPixel *, long, long, long, double, double,
                                    double);

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Requires buffering
  virtual bool RequiresBuffering();

public:

  /// Default constructor
  irtkEuclideanDistanceTransform(irtkDistanceTransformMode = irtkDistanceTransform3D);

  /// Destructor (empty).
  ~irtkEuclideanDistanceTransform() {};

  // Run distance transform
  virtual void Run();

  // Get Radial
  virtual void Radial();

  // Get Radial+Thickness
  virtual void TRadial();
};

template <class VoxelType> inline const char *irtkEuclideanDistanceTransform<VoxelType>::NameOfClass()
{
  return "irtkEuclideanDistanceTransform";
}

template <class VoxelType> inline bool irtkEuclideanDistanceTransform<VoxelType>::RequiresBuffering()
{
  return false;
}

#endif
