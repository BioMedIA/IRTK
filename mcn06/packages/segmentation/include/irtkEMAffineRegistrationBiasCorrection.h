/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKEMAFFINEREGISTRATIONBIASCORRECTION_H

#define _IRTKEMAFFINEREGISTRATIONBIASCORRECTION_H

#include <irtkImageFunction.h>
#include <irtkRegistration.h>
#include <irtkTransformation.h>
#include <irtkMLSimilarityMetric.h>
#include <irtkEMClassificationTemplateBiasCorrection.h>

/*

EM registration algorithm with bias field estimation and intensity matching

*/

class irtkEMAffineRegistrationBiasCorrection : public irtkEMClassificationTemplateBiasCorrection
{

  /// Uncorrected target image
  irtkGreyImage _reg_utr;

  /// Target image
  irtkGreyImage _reg_tr;

  /// Reference image
  irtkGreyImage _reg_rf;

  /// Reference image
  irtkRealImage _orig_rf;


  /// Registration
  irtkImageRegistration *_registration;

  ///Transformation
  irtkTransformation *_transformation;

  ///Interpolator
  irtkImageFunction *_interpolator;

  ///Image transformation
  irtkImageTransformation *_imagetransformation;

  irtkMLSimilarityMetric *_metric;


public:

  /// Estimates intensity matching
  virtual void RStep();

  /// Estimates bias field
  virtual void BStep();


public:

  ///Constructor
  irtkEMAffineRegistrationBiasCorrection(irtkRealImage &image, irtkRealImage &reference, int cP, irtkRealPixel padding, double voxelsize, char* parameters_name, irtkGreyImage& tr, irtkGreyImage& rf);
  /// Destructor
  ~irtkEMAffineRegistrationBiasCorrection();

  /// Intialize segmentation process
  virtual void Initialise();

  /// Execute one iteration and return log likelihood for GMM
  virtual double IterateGMM(int iteration);

  /// Bias correct target image
  virtual void CorrectTarget();

  /// Transform reference image
  void UpdateReference();

  void SaveTransformation(char*);

};

inline void irtkEMAffineRegistrationBiasCorrection::SaveTransformation(char* name)
{
  _transformation->Write(name);
}


#endif
