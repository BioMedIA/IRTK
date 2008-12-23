/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkEMClassificationBiasCorrection.h>

irtkEMClassificationBiasCorrection::irtkEMClassificationBiasCorrection() : irtkEMClassification()
{
  cerr<<"irtkEMClassificationBiasCorrection() ";
}

irtkEMClassificationBiasCorrection::irtkEMClassificationBiasCorrection(int noTissues, irtkRealImage **atlas) : irtkEMClassification(noTissues, atlas, NULL)
{
}

irtkEMClassificationBiasCorrection::irtkEMClassificationBiasCorrection(int noTissues, irtkRealImage **atlas, irtkRealImage *background) : irtkEMClassification(noTissues, atlas, background)
{
}

void irtkEMClassificationBiasCorrection::SetInput(const irtkRealImage &image)
{
  irtkEMClassification::SetInput(image);
  _uncorrected = image;
}

void irtkEMClassificationBiasCorrection::SetBiasField(irtkBiasField *biasfield)
{
  _biasfield = biasfield;
}

void irtkEMClassificationBiasCorrection::BStep()
{
  // Create bias correction filter
  _biascorrection.SetInput(&_uncorrected, &_estimate);
  _biascorrection.SetWeights(&_weights);
  _uncorrected.Write("uncorrected.nii.gz");
  _estimate.Write("estimate.nii.gz");
  _biascorrection.SetOutput(_biasfield);
  _biascorrection.SetPadding((short int) _padding);
  _biascorrection.Run();

  // Generate bias corrected image for next iteration
  _input = _uncorrected;
  _biascorrection.Apply(_input);
}

double irtkEMClassificationBiasCorrection::Iterate(int iteration)
{
  this->MStep();
  this->BStep();
  this->EStep();
  Print();
  return LogLikelihood();
}

double irtkEMClassificationBiasCorrection::IterateGMM(int iteration)
{
  if (iteration > 1) this->EStepGMM();
  this->MStepGMM();
  PrintGMM();
  this->WStep();
  this->BStep();
  _biasfield->Print();

  return LogLikelihoodGMM();
}

void irtkEMClassificationBiasCorrection::ApplyBias(irtkRealImage &image)
{
  _biascorrection.ApplyToImage(image);
}

void irtkEMClassificationBiasCorrection::ApplyBias(irtkGreyImage &image)
{
  _biascorrection.ApplyToImage(image);
}

void irtkEMClassificationBiasCorrection::ConstructBiasCorrectedImage(irtkRealImage &image)
{
  image = _input;
  _biascorrection.Apply(_input);
}
