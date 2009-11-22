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

#include <irtkEMAffineRegistrationBiasCorrection.h>

irtkEMAffineRegistrationBiasCorrection::irtkEMAffineRegistrationBiasCorrection(irtkRealImage &target, irtkRealImage &reference, int cP, irtkRealPixel padding, double voxelsize, char* parameters_name, irtkGreyImage& tr, irtkGreyImage& rf) : irtkEMClassificationTemplateBiasCorrection(target, reference, cP, padding, voxelsize)
{
  _reg_utr = tr;
  _reg_tr = tr;
  _reg_rf = rf;
  _orig_rf = reference;
  _registration = new irtkImageAffineRegistration;
  _transformation = new irtkAffineTransformation;
  //_registration = new irtkImageRigidRegistration;
  //_transformation = new irtkRigidTransformation;
  _registration->SetInput(&tr, &rf);
  _registration->SetOutput(_transformation);
  _registration->Read(parameters_name);
// _registration->SetTargetPadding((int) padding);

  _interpolator = new irtkLinearInterpolateImageFunction;
  _imagetransformation = new irtkImageTransformation;
  _imagetransformation->PutInterpolator(_interpolator);
  _imagetransformation->PutSourcePaddingValue(padding-1);
  _imagetransformation->PutTargetPaddingValue(padding-1);


}

irtkEMAffineRegistrationBiasCorrection::~irtkEMAffineRegistrationBiasCorrection()
{
  delete _registration;
  delete _transformation;
  delete _interpolator;
  delete _imagetransformation;
}

void irtkEMAffineRegistrationBiasCorrection::Initialise()
{
  irtkEMClassificationTemplateBiasCorrection::Initialise();
  _registration->classification=this;
  RStep();
}

void irtkEMAffineRegistrationBiasCorrection::RStep()
{
  GInit();
  _registration->Run();
  UpdateReference();
  Update();
}

void irtkEMAffineRegistrationBiasCorrection::BStep()
{
  irtkEMClassificationBiasCorrection::BStep();
  CorrectTarget();
  Update();
  _target.Write("_target.nii.gz");
}


double irtkEMAffineRegistrationBiasCorrection::IterateGMM(int iteration)
{
  LogLikelihoodGMM();
  if (iteration > 1) this->RStep();
  if (iteration > 1) this->IStep();
  if (iteration > 1) this->EStepGMM();
  this->MStepGMM();
  PrintGMM();
  this->WStep();
  this->BStep();
  _biasfield->Print();

  return LogLikelihoodGMM();
}

void irtkEMAffineRegistrationBiasCorrection::CorrectTarget()
{
  irtkEMClassificationTemplateBiasCorrection::CorrectTarget();
  ApplyBias(_reg_tr);
}

void irtkEMAffineRegistrationBiasCorrection::UpdateReference()
{
  _imagetransformation->SetInput(&_orig_rf, _transformation);
  irtkRealImage r(_orig_rf);
  r.Write("r.nii.gz");
  _imagetransformation->SetOutput(&r);
  cerr << _imagetransformation->GetTargetPaddingValue()<<endl;
  _imagetransformation->Run();
  cerr << _imagetransformation->GetTargetPaddingValue()<<endl;

  r.Write("regref.nii.gz");
  _reference=r;
}


