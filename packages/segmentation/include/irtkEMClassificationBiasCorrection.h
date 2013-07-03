/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKEMCLASSIFICATIONBIASCORRECTION_H

#define _IRTKEMCLASSIFICATIONBIASCORRECTION_H

#include <irtkImage.h>

#include <irtkProbabilisticAtlas.h>

#include <irtkEMClassification.h>

#include <irtkGaussianBlurring.h>

/*

EM maximisation algorithm with bias field estimation using Gaussian blurring and equal variance

*/

class irtkEMClassificationBiasCorrection : public irtkEMClassification
{

protected:
  /// Uncorrected image
  irtkRealImage _uncorrected;

  /// Bias field
  irtkRealImage _bias;
  
  /// Gaussian filter
  irtkGaussianBlurring<irtkRealPixel>* _gb;


 
public:

  /// Estimates bias field
  virtual void BStep();

  /// Estimates bias field
  virtual void WStep();

public:

  /// Constructor
  irtkEMClassificationBiasCorrection(int noTissues, irtkRealImage **atlas, double sigma);

  /// Constructor
  irtkEMClassificationBiasCorrection(double sigma);

  /// Constructor
  irtkEMClassificationBiasCorrection(int noTissues, irtkRealImage **atlas, irtkRealImage *background, double sigma);

  /// Destructor
  ~irtkEMClassificationBiasCorrection() ;

  /// Set image
  virtual void SetInput(const irtkRealImage &);

   /// Execute one iteration and return log likelihood
  virtual double Iterate(int iteration);

  /// Execute one iteration and return log likelihood for GMM
  virtual double IterateGMM(int iteration,bool equal_var, bool uniform_prior);

  /// Apply the bias 
  virtual void ApplyBias();

  /// Apply the bias to an image 
  virtual void ApplyBias(irtkRealImage &);

  /// Write resulting bias field
  void WriteBias(char* biasname);

};

inline void irtkEMClassificationBiasCorrection::WriteBias(char* biasname)
{
  _bias.Write(biasname);
}

#endif
