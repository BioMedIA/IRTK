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

#include <irtkMultiChannelImage.h>

irtkEMClassificationBiasCorrection::irtkEMClassificationBiasCorrection(double sigma) : irtkEMClassification()
{
  _gb = new irtkGaussianBlurring<irtkRealPixel>(sigma);
  //cerr<<"irtkEMClassificationBiasCorrection() ";
}

irtkEMClassificationBiasCorrection::irtkEMClassificationBiasCorrection(int noTissues, irtkRealImage **atlas, double sigma) : irtkEMClassification(noTissues, atlas)
{
 _gb = new irtkGaussianBlurring<irtkRealPixel>(sigma);
}

irtkEMClassificationBiasCorrection::irtkEMClassificationBiasCorrection(int noTissues, irtkRealImage **atlas, irtkRealImage *background, double sigma) : irtkEMClassification(noTissues, atlas, background)
{
 _gb = new irtkGaussianBlurring<irtkRealPixel>(sigma);
}

irtkEMClassificationBiasCorrection::~irtkEMClassificationBiasCorrection()
{
  //irtkEMClassification::~irtkEMClassification();
  delete _gb;
}

void irtkEMClassificationBiasCorrection::SetInput(const irtkRealImage &image)
{
  irtkEMClassification::SetInput(image);
  _uncorrected = image;
  _bias=image;
  _bias.PutMinMax(0,0);
}

void irtkEMClassificationBiasCorrection::BStep()
{
  int i;
  double scale = 1000;
  irtkRealImage residual(_input);
  irtkRealImage wresidual(_input);
  // Because of equal sigmas mask is just normalized version of weights

  cerr<<"Calculating bias ...";

  //calculate residual image
  irtkRealPixel *pi=_input.GetPointerToVoxels();
  irtkRealPixel *pw=_weights.GetPointerToVoxels();
  irtkRealPixel *pe=_estimate.GetPointerToVoxels();
  irtkRealPixel *pm=_mask.GetPointerToVoxels();
  irtkRealPixel *pr=residual.GetPointerToVoxels();
  irtkRealPixel *prw=wresidual.GetPointerToVoxels();

  //_output.First();
  //_atlas.First();

  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    if ((*pm == 1)&&(*pi != _padding)) {
      *pr=*pi / *pe;
      *prw= (*pw) * log(*pr)*scale;
    } else {
      *pr=_padding;
      *prw=_padding;
    }
    pi++;
    pm++;
    pw++;
    pr++;
    prw++;
    pe++;
    //_output.Next();
    //_atlas.Next();
  }
  residual.Write("residual.nii.gz");
  wresidual.Write("wresidual.nii.gz");
  _weights.Write("_weights.nii.gz");
  _bias.Write("bias-start.nii.gz");
  
  irtkMultiChannelImage mch;
  mch.AddImage(wresidual);
  mch.AddImage(_weights);
  mch.SetMask(_mask);
  mch.SetPadding(0);
  //mch.Log(0);
  //mch.Log(1);
  //mch.Write(0,"logresidual.nii.gz");
  //mch.Write(1,"logweights.nii.gz");

  _gb->SetInput(&mch.GetImage(0));
  _gb->SetOutput(&mch.GetImage(0));
  _gb->Run();
  mch.Write(0,"wresidualblurred.nii.gz");
  
  _gb->SetInput(&mch.GetImage(1));
  _gb->SetOutput(&mch.GetImage(1));
  _gb->Run();
  mch.Write(1,"weights-blurred.nii.gz");

  //calculate weighted blurring of log residual
  irtkRealImage res;
  res = mch.Divide();
  mch.SetImage(0,res);
  mch.Write(0,"logresidualblurredweigted.nii.gz");

  //Calculate bias
  //mch.Brainmask();
  _bias += mch.GetImage(0);
  //_bias = mch.GetImage(0);
  _bias.Write("bias.nii.gz");
  mch.Write(0,"diffbias.nii.gz");
  
  //set the mean of the bias field to zero
  irtkRealPixel *pb=_bias.GetPointerToVoxels();
  pi=_input.GetPointerToVoxels();
  pm=_mask.GetPointerToVoxels();
  double sum=0;
  int num=0;
  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    if ((*pm == 1)&&(*pi != _padding)) {
      sum+=*pb;
      num++;
    } 
    pi++;
    pm++;
    pb++;
  }
  double mean = sum/num;
  //irtkRealImage diffbias=mch.GetImage(0);

  irtkRealPixel *pd=mch.GetImage(0).GetPointerToVoxels();
  pb=_bias.GetPointerToVoxels();
  pi=_input.GetPointerToVoxels();
  pm=_mask.GetPointerToVoxels();
  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    //if ((*pm == 1)&&(*pi != _padding)) {
      *pb-=mean;
      *pd-=mean;
    //} 
    pi++;
    pm++;
    pb++;
    pd++;
  }
  cerr<<"Adjusted mean of the bias "<<mean<<" to zero."<<endl;
  _bias.Write("zerobias.nii.gz");
  
  
  mch.Exp(0,1000);  
  mch.Write(0,"residualblurredweigted.nii.gz");
  mch.Brainmask();
  mch.Write(0,"expbias.nii.gz");

  //Correct input
  mch.SetImage(1,mch.GetImage(0));
  mch.SetImage(0,_input);
  //mch.Exp(1,1);
  mch.Write(0,"input.nii.gz");
  //mch.Write(1,"biasexp.nii.gz");
  _input=mch.Divide();
  _input.Write("corrected.nii.gz");
  cerr<<"done."<<endl;
}


double irtkEMClassificationBiasCorrection::Iterate(int iteration)
{
  this->EStep();
  this->MStep();
  this->WStep();
  this->BStep();

  Print();
  return LogLikelihood();
}

double irtkEMClassificationBiasCorrection::IterateGMM(int iteration, bool equal_var, bool uniform_prior)
{
  if (iteration > 1) this->EStepGMM();
  //this->MStepGMM();
  if (equal_var) this->MStepVarGMM(uniform_prior);
  else this->MStepGMM(uniform_prior);
  PrintGMM();
  this->WStep();
  this->BStep();
  PrintGMM();

  return LogLikelihoodGMM();
}

void irtkEMClassificationBiasCorrection::WStep()
{
  int i,k;
  double num, den;
  cerr<<"Calculating weights ...";
  irtkRealPixel *pi=_input.GetPointerToVoxels();
  irtkRealPixel *pw=_weights.GetPointerToVoxels();
  irtkRealPixel *pwB=_weightsB.GetPointerToVoxels();
  irtkRealPixel *pb=_bias.GetPointerToVoxels();
  irtkRealPixel *pe=_estimate.GetPointerToVoxels();
  irtkRealPixel *pm=_mask.GetPointerToVoxels();
  _output.First();
  _atlas.First();

  for (i=0; i< _input.GetNumberOfVoxels(); i++) {
    if ((*pm == 1)&&(*pi != _padding)) {
      num=0;
      den=0;
      for (k=0; k<_number_of_tissues; k++) {
        num += _output.GetValue(k)*_mi[k]/_sigma[k];
        den += _output.GetValue(k)/_sigma[k];
      }
      *pw=den*(*pi);
      *pwB=_output.GetValue(0)/_sigma[0];
      *pe=num/den;

      /*if(_background_tissue >=0)
      {
        if(_atlas.GetValue(_background_tissue) == 1) *pe=_padding;
      }
      */
    } else {
      *pw=_padding;
      *pe=_padding;
    }

    pi++;
    pm++;
    pw++;
    pb++;
    pwB++;
    pe++;
    _output.Next();
    _atlas.Next();
  }
  _estimate.Write("estimate.nii.gz");
  //_weights.Write("_weights.nii.gz");
  //_weightsR.Write("_weightsR.nii.gz");
  //_weightsB.Write("_weightsB.nii.gz");
  cerr<<"done."<<endl;
}


void irtkEMClassificationBiasCorrection::ApplyBias()
{
  
}

void irtkEMClassificationBiasCorrection::ApplyBias(irtkRealImage &image)
{
  irtkMultiChannelImage mch;
  mch.AddImage(image);
  mch.SetPadding(0);
  mch.CreateMask();
  mch.AddImage(_bias);
  mch.Exp(1,1000);
  mch.Brainmask();
  image = mch.Divide();
 
}

