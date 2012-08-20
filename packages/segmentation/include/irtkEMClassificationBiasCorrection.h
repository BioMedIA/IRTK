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

#include <irtkBiasField.h>

#include <irtkBiasCorrection.h>

/*

EM maximisation algorithm with bias field estimation

*/

class irtkEMClassificationBiasCorrection : public irtkEMClassification
{

protected:
  /// Uncorrected image
  irtkRealImage _uncorrected;

  /// Bias field
  irtkBiasField *_biasfield;

  /// Bias field correction filter
  irtkBiasCorrection _biascorrection;

public:

  /// Estimates bias field
  virtual void BStep();

public:

  /// Constructor
  irtkEMClassificationBiasCorrection(int noTissues, irtkRealImage **atlas);

  /// Constructor
  irtkEMClassificationBiasCorrection();

  /// Constructor
  irtkEMClassificationBiasCorrection(int noTissues, irtkRealImage **atlas, irtkRealImage *background);

  /// Set image
  virtual void SetInput(const irtkRealImage &);

  /// Set bias field
  virtual void SetBiasField(irtkBiasField *);

  /// Execute one iteration and return log likelihood
  virtual double Iterate(int iteration);

  /// Execute one iteration and return log likelihood for GMM
  virtual double IterateGMM(int iteration, bool equal_var, bool uniform_prior);

  /// Compute the bias corrected image
  virtual void ConstructBiasCorrectedImage(irtkRealImage &);

  /// Apply the bias to an image
  virtual void ApplyBias(irtkRealImage &);

  /// Apply the bias to Grey image icluding logarithmic transform
  virtual void ApplyBias(irtkGreyImage &);

  /// Write resulting bias field
  void WriteBias(char* biasname);

};

inline void irtkEMClassificationBiasCorrection::WriteBias(char* biasname)
{
  _biasfield->Write(biasname);
}

#endif
