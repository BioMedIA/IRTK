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

#include <irtkEMClassificationTemplateBiasCorrection.h>
#include <irtkHistogram.h>
#include <irtkMultiChannelImage.h>
#include <irtkResampling.h>
#include <irtkResamplingWithPadding.h>


irtkEMClassificationTemplateBiasCorrection::irtkEMClassificationTemplateBiasCorrection(irtkRealImage &target, irtkRealImage &reference, double spacing, irtkRealPixel padding, double voxelsize)
{
  cerr<<endl;
  cerr<<"irtkEMClassificationTemplateBiasCorrection: "<<endl;
  cerr<<"spacing="<<spacing<<" padding="<<padding<<" voxelsize="<<voxelsize<<endl;

  _uncorrected_target = target;
  _target=target;
  _reference=reference;

  _padding=padding;
  _voxelsize=voxelsize;
  _number_of_voxels = 0;
  _init=false;


  //IntensityInit();
  irtkMultiChannelImage mch;
  mch.SetPadding((int) _padding);
  mch.AddImage(Resample(_uncorrected_target));
  mch.AddImage(Resample(_reference));
  mch.Log(0);
  mch.Log(1);
  mch.CreateMask();
  mch.Brainmask();


  _d_uncorrected_target = mch.GetImage(0);
  _d_target = _d_uncorrected_target;
  _d_reference = mch.GetImage(1);
  _d_rm = _d_reference;

  SetInput(mch.Subtract());

  _biasfield = new irtkBSplineBiasField(_d_target, spacing, spacing, spacing);
  PutBoundingBox();
  _biascorrection.SetOutput(_biasfield);
   //CorrectTarget();
   //Update();


}

irtkEMClassificationTemplateBiasCorrection::~irtkEMClassificationTemplateBiasCorrection()
{
  delete _biasfield;
}

void irtkEMClassificationTemplateBiasCorrection::PutBoundingBox()
{

  double x1,y1,z1,x2,y2,z2;
  int treshold = 0;
  // Default roi
  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = _input.GetX();
  y2 = _input.GetY();
  z2 = _input.GetZ();

  int i,j,k;
/////////z coordinate
  int sum=0;
  for(k=_input.GetZ()-1;k>=0;k--)
  {
    sum=0;
    for(j=_input.GetY()-1;j>=0;j--)
      for(i=_input.GetX()-1;i>=0;i--)
      {
        if(_mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  z2=k;
  cerr<<z2<<endl;

  sum=0;
  for(k=0; k<=_input.GetZ()-1;k++)
  {
    sum=0;
    for(j=_input.GetY()-1;j>=0;j--)
      for(i=_input.GetX()-1;i>=0;i--)
      {
        if(_mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  z1=k;
  cerr<<z1<<endl;

/////////y coordinate
  sum=0;
  for(j=_input.GetY()-1;j>=0;j--)
  {
    sum=0;
    for(k=_input.GetZ()-1;k>=0;k--)
      for(i=_input.GetX()-1;i>=0;i--)
      {
        if(_mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  y2=j;
  cerr<<y2<<endl;

  sum=0;
  for(j=0; j<=_input.GetY()-1;j++)
  {
    sum=0;
    for(k=_input.GetZ()-1;k>=0;k--)
      for(i=_input.GetX()-1;i>=0;i--)
      {
        if(_mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  y1=j;
  cerr<<y1<<endl;

/////////x coordinate
  sum=0;
  for(i=_input.GetX()-1;i>=0;i--)
  {
    sum=0;
    for(k=_input.GetZ()-1;k>=0;k--)
      for(j=_input.GetY()-1;j>=0;j--)
      {
        if(_mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  x2=i;
  cerr<<x2<<endl;

  sum=0;
  for(i=0; i<=_input.GetX()-1;i++)
  {
    sum=0;
    for(k=_input.GetZ()-1;k>=0;k--)
      for(j=_input.GetY()-1;j>=0;j--)
      {
        if(_mask.Get(i,j,k)>0) sum++;
      }
    if (sum>treshold) break;
  }

  x1=i;
  cerr<<x1<<endl;

  cerr<<x1<<" "<<y1<<" "<<z1<<" "<<x2<<" "<<y2<<" "<<z2<<endl;
  x1 -= 1; 
  y1 -= 1;
  z1 -= 1;
  x2 += 1; 
  y2 += 1;
  z2 += 1;
  
  _input.ImageToWorld(x1, y1, z1);
  _input.ImageToWorld(x2, y2, z2);
  _biasfield->PutBoundingBox(x1,y1,z1,x2,y2,z2);


}

void irtkEMClassificationTemplateBiasCorrection::Initialise(bool nomatch)
{
  if (!nomatch) IStep();
  //IntensityInit();
  SetInput(_uncorrected);
  InitialiseGMMParameters();
  PrintGMM();
  InitialisePosteriors();
  LogLikelihoodGMM();
}

void irtkEMClassificationTemplateBiasCorrection::InitialiseGMMParameters()
{
  int i;

  cerr <<"Estimating GMM parameters ... ";
  irtkHistogram h(100);
  irtkRealPixel imin, imax;
  _input.GetMinMaxPad(&imin, &imax,_padding);
  cerr<<" min = "<<imin<<", max = "<<imax<<" ... ";
  h.PutMin(imin);
  h.PutMax(imax);
  irtkRealPixel *ptr=_input.GetPointerToVoxels();
  for (i=0; i<_input.GetNumberOfVoxels(); i++) {
    if (*ptr!=_padding) h.AddSample(*ptr);
    ptr++;
  }
  double mean, variance;
  //mean = h.Mean();
  mean = 0;
  variance=h.Variance();
  cerr<<"mean="<<mean<<" variance="<<sqrt(variance)<<" ... done."<<endl;

  _number_of_tissues=2;
  _mi=new double[2];
  _sigma=new double[2];
  _c=new double[2];
  _mi[0]=_mi[1]=mean;
  _sigma[0]=variance/4;
  _sigma[1]=variance*4;
  _c[0]=_c[1]=0.5;
}


void irtkEMClassificationTemplateBiasCorrection::IStep()
{
  irtkRealPixel *pt, *pr, *pw=NULL, *pm;
  int i,n=0;
  double *x,*y,*w;
  irtkRealPixel xmin, xmax, ymin, ymax;

  //initialise
  _d_target.GetMinMax(&ymax,&ymin);
  _d_reference.GetMinMax(&xmax,&xmin);

  pt = _d_target.GetPointerToVoxels();

  for ( i=0; i<_d_target.GetNumberOfVoxels(); i++) {
    if (*pt != _padding) n++;
    pt++;
  }

  x = new double[n];
  y = new double[n];
  w = new double[n];

  pt = _d_target.GetPointerToVoxels();
  pr = _d_reference.GetPointerToVoxels();
  if (_init) pw = _weights.GetPointerToVoxels();

  n=0;
  for ( i=0; i<_d_target.GetNumberOfVoxels(); i++) {
    if (*pt != _padding) {
      x[n] = *pr;
      y[n] = *pt;
      if (_init) w[n]=*pw;
      else w[n]=1;
      n++;
      if (*pr<xmin) xmin=*pr;
      if (*pr>xmax) xmax=*pr;
      if (*pt<ymin) ymin=*pt;
      if (*pt>ymax) ymax=*pt;
    }
    pt++;
    pr++;
    if (_init) pw++;
  }

  _matching.SetMinMax(xmin,xmax);
  _matching.SetYMinMax(ymin,ymax);
  //_matching.WeightedPCA(x,y,w,n);
  //_matching.WeightedLeastSqaures(x,y,w,n);
  _matching.MatchMeanAndVariance(x,y,w,n);

  pr = _d_reference.GetPointerToVoxels();
  pm = _d_rm.GetPointerToVoxels();

  for ( i=0; i<_d_target.GetNumberOfVoxels(); i++) {
    if (*pr != _padding) *pm=_matching.Lin(*pr);
    if (*pm < 0) *pm = 0;
    pr++;
    pm++;
  }
  _d_rm.Write("rm.nii.gz");


  irtkMultiChannelImage mch;
  mch.SetPadding((int)_padding);
  mch.AddImage(_d_uncorrected_target);
  mch.AddImage(_d_rm);

  _uncorrected=mch.Subtract();

  mch.SetImage(0,_d_target);
  _input=mch.Subtract();


  delete[] x;
  delete[] y;
  delete[] w;
}

void irtkEMClassificationTemplateBiasCorrection::Update()
{
  irtkMultiChannelImage mch;
  mch.SetPadding((int) _padding);
  mch.AddImage(_uncorrected_target);
  mch.AddImage(_reference);
  mch.CreateMask();
  mch.Brainmask();
  mch.Log(0);
  mch.Log(1);

  _d_uncorrected_target = Resample(mch.GetImage(0));
  _d_uncorrected_target.Write("_dut.nii.gz");
  _d_target = _d_uncorrected_target;
  _d_target.Write("_dt.nii.gz");
  ApplyBias(_d_target);
  _d_target.Write("_dt.nii.gz");

  _d_reference = Resample(mch.GetImage(1));
  _d_reference.Write("_dr.nii.gz");
  _d_rm = _d_reference;
  //MatchIntensity(_d_rm);
  _d_rm.Write("_drm.nii.gz");


}

irtkRealImage irtkEMClassificationTemplateBiasCorrection::Resample( irtkRealImage& image)
{
  irtkRealImage i3;
  cerr << "Resampling with padding ... ";
  irtkResamplingWithPadding<irtkRealPixel> resampling(_voxelsize, _voxelsize, _voxelsize, _padding);
  resampling.SetInput(&image);
  resampling.SetOutput(&i3);
  resampling.Run();
  cerr<<"done."<<endl;
  return i3;
}

double irtkEMClassificationTemplateBiasCorrection::IterateGMM(int iteration)
{
  //if (iteration > 1) this->IStep();
  if (iteration > 1) this->EStepGMM();
  this->MStepGMM();
  _mi[0]=0;
  _mi[1]=0;
  PrintGMM();
  this->WStep();
  this->BStep();
  //Update();
  _biasfield->Print();

  return LogLikelihoodGMM();
}

void irtkEMClassificationTemplateBiasCorrection::CorrectTarget()
{
  irtkMultiChannelImage corrected;
  corrected.SetPadding(0);
  corrected.AddImage(_uncorrected_target);
  corrected.AddImage(_uncorrected_target);
  corrected.Log(0);
  ApplyBias(corrected.GetImage(0));
  corrected.Exp(0);
  //corrected.AdjustMean(0,1);
  _target = corrected.GetImage(0);
}


double irtkEMClassificationTemplateBiasCorrection::MatchIntensity(double x)
{
  return _matching.Lin(x);
}

void irtkEMClassificationTemplateBiasCorrection::MatchIntensity(irtkRealImage &image)
{
  irtkRealPixel *p=image.GetPointerToVoxels();
  for ( int i=0; i<image.GetNumberOfVoxels(); i++) {
    if (*p != _padding) *p=_matching.Lin(*p);
    if (*p < 0) *p = 0;
    p++;
  }
}

double irtkEMClassificationTemplateBiasCorrection::PointLogLikelihoodGMM(double x, double y)
{

  double lx, ly;

  if (x<0) return 0;
  if (x==0) x=1;
  if (y<=0) y=1;

  lx=1000*log(x);
  ly=MatchIntensity(1000*log(y));
  if (ly<0) ly=0;
  return irtkEMClassification::PointLogLikelihoodGMM(lx-ly);
}

double irtkEMClassificationTemplateBiasCorrection::PointLogLikelihoodGMMnomatch(double x, double y)
{
  return irtkEMClassification::PointLogLikelihoodGMM(x-y);
}


void irtkEMClassificationTemplateBiasCorrection::IntensityInit()
{

  cerr<<"IntensityInit: "<<endl;
  irtkRealPixel *pt, *pr;
  int i,n=0;
  double *x,*y,*w;
  irtkRealPixel xmin, xmax, ymin, ymax;

  //initialise
  _target.GetMinMax(&ymax,&ymin);
  _reference.GetMinMax(&xmax,&xmin);

  
  pt = _target.GetPointerToVoxels();
  pr = _reference.GetPointerToVoxels();
  _target.Write("target.nii.gz");
  _reference.Write("reference.nii.gz");
  
  for( i=0; i<_target.GetNumberOfVoxels(); i++) 
  {
    if((*pt != _padding)&&(*pr != _padding)) n++;
    pt++;
    pr++;
  }
  
  x = new double[n];
  y = new double[n];
  w = new double[n];

  pt = _target.GetPointerToVoxels();
  pr = _reference.GetPointerToVoxels();
  //if (_init) pw = _weights.GetPointerToVoxels();

  n=0;
  for( i=0; i<_target.GetNumberOfVoxels(); i++) 
  {
    if((*pt != _padding)&&(*pr != _padding))
    {
      x[n] = *pr;
      y[n] = *pt;
      //if (_init) w[n]=*pw;
      w[n]=1;
      n++;
      if(*pr<xmin) xmin=*pr;
      if(*pr>xmax) xmax=*pr;
      if(*pt<ymin) ymin=*pt;
      if(*pt>ymax) ymax=*pt;
    }
    pt++;
    pr++;
  }

  _matching.SetMinMax(xmin,xmax);
  _matching.SetYMinMax(ymin,ymax);
  
  _matching.MatchMeanAndVariance(x,y,w,n);
  
   cerr<<endl<<"Updating reference ... "<<endl;
    pr = _reference.GetPointerToVoxels();
    for( i=0; i<_reference.GetNumberOfVoxels(); i++) 
    {
      if(*pr != _padding)
      {
        *pr=_matching.Lin(*pr);
      }
      pr++;
   }
   _reference.Write("matched-reference.nii.gz");

  delete[] x;
  delete[] y;
  delete[] w;
}

