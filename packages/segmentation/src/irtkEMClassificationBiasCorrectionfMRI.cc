/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkEMClassificationBiasCorrection.cc 933 2013-07-03 15:33:19Z mm3 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2013-07-03 16:33:19 +0100 (Wed, 03 Jul 2013) $
  Version   : $Revision: 933 $
  Changes   : $Author: mm3 $

=========================================================================*/

#include <irtkImage.h>

#include <irtkEMClassificationBiasCorrectionfMRI.h>

#include <irtkMultiChannelImage.h>


/*
void irtkEMClassificationBiasCorrectionfMRI::SetInput(const irtkRealImage &image)
{
  irtkEMClassification::SetInput(image);
  _uncorrected = image;
  
  irtkImageAttributes attr = image.GetImageAttributes();
  attr._t = 1;
  irtkRealImage bias(attr);
  _bias=bias;
  _bias=0;
}
*/

irtkEMClassificationBiasCorrectionfMRI::irtkEMClassificationBiasCorrectionfMRI(double sigma) : 
  irtkEMClassificationBiasCorrection(sigma)
{
}


void irtkEMClassificationBiasCorrectionfMRI::BStep()
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
  //residual.Write("residual.nii.gz");
  //wresidual.Write("wresidual.nii.gz");
  //_weights.Write("_weights.nii.gz");
  //_bias.Write("bias-start.nii.gz");
  
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
  //mch.Write(0,"wresidualblurred.nii.gz");
  
  _gb->SetInput(&mch.GetImage(1));
  _gb->SetOutput(&mch.GetImage(1));
  _gb->Run();
  //mch.Write(1,"weights-blurred.nii.gz");

  //calculate weighted blurring of log residual
  irtkRealImage res;
  res = mch.Divide();
  mch.SetImage(0,res);
  //mch.Write(0,"logresidualblurredweigted.nii.gz");
  NormalizeInTime(mch.GetImage(0));

  //Calculate bias
  //mch.Brainmask();
  _bias += mch.GetImage(0);
  //_bias = mch.GetImage(0);
  //_bias.Write("bias.nii.gz");
  //mch.Write(0,"diffbias.nii.gz");
  
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
  _bias.Write("bias.nii.gz");
  
  
  mch.Exp(0,1000);  
  //mch.Write(0,"residualblurredweigted.nii.gz");
  mch.Brainmask();
  //mch.Write(0,"expbias.nii.gz");

  //Correct input
  mch.SetImage(1,mch.GetImage(0));
  mch.SetImage(0,_input);
  //mch.Exp(1,1);
  //mch.Write(0,"input.nii.gz");
  //mch.Write(1,"biasexp.nii.gz");
  _input=mch.Divide();
  //_input.Write("corrected.nii.gz");
  cerr<<"done."<<endl;
}

void irtkEMClassificationBiasCorrectionfMRI::NormalizeInTime(irtkRealImage &image)
{
  int i,j,k,t;
  double num,sum,average;
  
  //image.Write("nit-image.nii.gz");
  for(i=0; i<image.GetX(); i++)
    for(j=0; j<image.GetY(); j++)
      for(k=0; k<image.GetZ(); k++)
      {
	num=0;
	sum=0;
	for(t=0;t<image.GetT();t++)
	{
	  num+=image(i,j,k,t);
	  sum++;
	}
	average = num/sum;
	for(t=0;t<image.GetT();t++)
	  image(i,j,k,t)=average;
      }
  //image.Write("nit-imagenorm.nii.gz");
}