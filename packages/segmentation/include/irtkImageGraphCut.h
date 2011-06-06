/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _irtkImageGraphCut_H

#define _irtkImageGraphCut_H

#include <irtkImage.h>
#include <GCoptimization.h>
/**
 * Abstract base class for any general image to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image as output. Each
 * derived class has to implement all abstract member functions.
 */

template <class VoxelType> class irtkImageGraphCut : public irtkObject
{

protected:

  /// Number of input images
  int _numberOfImages;

  /// spacing
  double _dx,_dy,_dz,_dt;

  /// Number of Weight images
  double _numberOfWeights;

  /// Input image for filter
  irtkGenericImage<VoxelType> **_input;

  /// Output image for filter
  irtkGenericImage<VoxelType> *_output;

  /// source weight image
  irtkRealImage **_weight;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /** Finalize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Finalize();

  /// add regionweight to graph
  virtual void AddBoundaryTerm(GCoptimizationGeneralGraph* graph, int count, 
	  int i,int j, int k, int l,
	  int xoff, int yoff, int zoff, int toff, double divide);


  /// geometry mode 0 no connection 1 x 2 xy 3 xyz 4 xyzt
  int _mode;

  /// number of labels
  int _labels;

public:

  /// Constructor
  irtkImageGraphCut();

  /// Deconstuctor
  virtual ~irtkImageGraphCut();

  /// Set input image for filter
  virtual void SetInput (int, irtkGenericImage<VoxelType> **,int, irtkRealImage **);

  /// Set input image for filter
  virtual void SetInput (irtkGenericImage<VoxelType> *,int, irtkRealImage **);

  /// Set output image for filter
  virtual void SetOutput(irtkGenericImage<VoxelType> *);

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

template <class VoxelType> const char *irtkImageGraphCut<VoxelType>::NameOfClass()
{
  return "irtkImageGraphCut";
}

template <class VoxelType> void irtkImageGraphCut<VoxelType>::SetMode(int i)
{
  _mode = i;
}

template <class VoxelType> VoxelType irtkImageGraphCut<VoxelType>::MinValue(VoxelType first, VoxelType second)
{
  if(first > second)
	  return second;
  else
	  return first;
}

template <class VoxelType> VoxelType irtkImageGraphCut<VoxelType>::MaxValue(VoxelType first, VoxelType second)
{
  if(first < second)
	  return second;
  else
	  return first;
}

#endif
