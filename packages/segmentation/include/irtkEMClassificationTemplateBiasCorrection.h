/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKEMCLASSIFICATIONTAMPLATEBIASCORRECTION_H

#define _IRTKEMCLASSIFICATIONTAMPLATEBIASCORRECTION_H

#include <irtkImage.h>

#include <irtkEMClassificationBiasCorrection.h>

#include <irtkIntensityMatching.h>


/*

EM maximisation algorithm with bias field estimation and intensity matching

*/

class irtkEMClassificationTemplateBiasCorrection : public irtkEMClassificationBiasCorrection
{

protected:

  /// Uncorrected target image
  irtkRealImage _uncorrected_target;

  /// Target image
  irtkRealImage _target;

  /// Reference image
  irtkRealImage _reference;

  ///Downsampled uncorrected image
  irtkRealImage _d_uncorrected_target;

  /// Downsamped target image
  irtkRealImage _d_target;

  /// Downsampled reference image
  irtkRealImage _d_reference;

  /// Downsampled intensity-matched reference image
  irtkRealImage _d_rm;


  /// Intensity matching
  irtkIntensityMatching _matching;

  /// Voxel size for downsampled image
  double _voxelsize;
  /// Indicates whether initialization of all members is fully finished
  bool _init;

public:

  /// Estimates intensity matching
  virtual void IStep();

public:

  ///Constructor
  irtkEMClassificationTemplateBiasCorrection(irtkRealImage &image, irtkRealImage &reference, double spacing, irtkRealPixel padding, double voxelsize);
  /// Destructor
  ~irtkEMClassificationTemplateBiasCorrection();

  /// Intialize segmentation process
  void Initialise();

  /// Estimate initial GMM parameters
  void InitialiseGMMParameters();

  /// Returns resampled image
  irtkRealImage Resample( irtkRealImage& image);

  /// Execute one iteration and return log likelihood for GMM
  virtual double IterateGMM(int iteration);

  /// Bias correct target image
  void CorrectTarget();

  /// Write target image
  void WriteTarget( char*);

  /// Write reference image
  void WriteReference(char*);

  ///convert intensity in reference image to intensity in target image
  virtual double MatchIntensity(double x);

  virtual double PointLogLikelihoodGMM(double x, double y);
  virtual double PointLogLikelihoodGMMnomatch(double x, double y);

  virtual void Update();
  void MatchIntensity(irtkRealImage&);
  void PutBoundingBox();
  void IntensityInit();




};

inline void irtkEMClassificationTemplateBiasCorrection::WriteTarget(char *name)
{
  _target.Write(name);
}

inline void irtkEMClassificationTemplateBiasCorrection::WriteReference(char *name)
{
  _reference.Write(name);
}


#endif
