/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkDynamicSegmentation.h 2 2008-12-23 12:40:14Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2008-12-23 12:40:14 +0000 (‰∫? 23 ÂçÅ‰∫åÊú?2008) $
  Version   : $Revision: 2 $
  Changes   : $Author: dr $

=========================================================================*/

#ifndef _irtkDynamicSegmentation_H

#define _irtkDynamicSegmentation_H

/**
 * Abstract base class for any general image to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image as output. Each
 * derived class has to implement all abstract member functions.
 */

template <class VoxelType> class irtkDynamicSegmentation : public irtkObject
{

protected:

  /// Input image (intensity)
  irtkGenericImage<VoxelType> *_input;

  /// Multiple image pointer
  irtkGenericImage<VoxelType> **_multiptr;

  /// Output image (label)
  irtkGenericImage<VoxelType> *_output;

  /// Input transformation
  irtkMultiLevelFreeFormTransformation *_mffd;

  /// Number of Labels
  int _NumberofLabels;

  /// Mode (Multi or Single Graphcut)
  int _mode;

  irtkRealImage **_atlas;

  /// 3D Gaussian image
  irtkGaussian **_gaussian;

  /// Gaussian norm
  double **_normal;

  /// motion maginitude
  irtkGenericImage<VoxelType> *_motion;

  /// jacobian determine
  irtkGenericImage<VoxelType> *_jacobian;

  /// probability
  irtkProbabilisticAtlas _probability;

  /// Graphu cut optimizer
  irtkImageGraphCut<VoxelType> *_graphcut;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /** Finalize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Finalize();

public:

  /// Constructor
  irtkDynamicSegmentation();

  /// Deconstuctor
  virtual ~irtkDynamicSegmentation();

  /// Set input image for filter
  virtual void SetInput (irtkGenericImage<VoxelType> *,irtkMultiLevelFreeFormTransformation *,  irtkRealImage **, int);

  /// Set output image for filter
  virtual void SetOutput(irtkGenericImage<VoxelType> *);

  /// Set update
  virtual void UpdateMFFD(int = 1);

  /// Run filter on entire image lambda is the coefficent for region term.
  virtual void Run();

  /// Returns the name of the class
  virtual const char *NameOfClass();
};

template <class VoxelType> const char *irtkDynamicSegmentation<VoxelType>::NameOfClass()
{
  return "irtkDynamicSegmentation";
}

#endif
