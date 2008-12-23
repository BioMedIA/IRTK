/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkIntensityMatching.h>

#include <irtkImage.h>

irtkIntensityMatching::irtkIntensityMatching()
{
  _min=0;
  _max=0;
  _m=1;
  _b=0;
}

irtkIntensityMatching::irtkIntensityMatching(double min, double max, double m=1, double b=0)
{
  _min=min;
  _max=max;
  _m=m;
  _b=b;
}

irtkIntensityMatching::irtkIntensityMatching(const irtkIntensityMatching &f)
{
  _m=f._m;
  _b=f._b;
  _min=f._min;
  _max=f._max;
}


irtkIntensityMatching::~irtkIntensityMatching()
{
  _min=0;
  _max=0;
  _m=1;
  _b=0;
}



void irtkIntensityMatching::Print()
{
  // Write no. of control points
  cerr<<"y = "<<_m<<" * x + "<<_b<<endl;
  cout << "Interval: <" << _min << ", " << _max<<">"<< endl;
}



void irtkIntensityMatching::PCA(double *x1, double *y1, double *weights, int no)
{
  //Calculates the largest 2D component by PCA - as the mean + eigenvector belonging to bigger eigenvalue of the covariance matrix
  cerr<<"Start principal component regression"<<endl;

  irtkGreyImage hist(101,101,1), trans(101,101,1);
  double xxx,yyy;
  int i;

  for (i=0; i<no; i++) {
    xxx=100*(x1[i]-_min)/(_max-_min);
    yyy=100*(y1[i]-_ymin)/(_ymax-_ymin);
    if ((xxx>=0)&&(xxx<=100)&&(yyy>=0)&&(yyy<=100))hist.Put(round(xxx),round(yyy),0,hist.Get(round(xxx),round(yyy),0)+1);
  }
  _hist=hist;
  hist.Write("hist.nii.gz");

  //expectation values and covarinaces
  double x=0,y=0,xy=0,x2=0,y2=0;
  double ex,ey,ex2,exy,ey2;
  double varx, vary, covxy;
  double k,l;

  for (i=0; i<no; i++) {
    x += x1[i];
    y += y1[i];
    xy += x1[i]*y1[i];
    x2 += x1[i]*x1[i];
    y2 += y1[i]*y1[i];
  }

  ex  = x/no;
  ey  = y/no;
  ex2 = x2/no;
  exy = xy/no;
  ey2 = y2/no;

  varx  = ex2-ex*ex;
  vary  = ey2-ey*ey;
  covxy = exy-ex*ey;

  k = (varx-vary)/(2*covxy);
  l = k+sqrt(1+k*k);

  _m = 1/l;
  _b = ey-ex/l;

  cerr<<"y = "<<_m<<" * x + "<<_b<<endl;

  double origx;
  for (i=0; i<=100; i++) {
    origx=_min+(_max-_min)*(i/100.0);
    yyy=100*(Lin(origx)-_ymin)/(_ymax-_ymin);
    if ((yyy>=0)&&(yyy<=100))trans.Put(i,round(yyy),0,1);
  }
  trans.Write("trans.nii.gz");
  _trans=trans;




}

void irtkIntensityMatching::WeightedPCA(double *x1, double *y1, double *weights, int no)
{
  //Calculates the largest 2D component by PCA - as the mean + eigenvector belonging to bigger eigenvalue of the covariance matrix
  cerr<<"Start weighted principal component regression"<<endl;

  irtkGreyImage hist(101,101,1), trans(101,101,1);
  double xxx,yyy;
  int i;

  for (i=0; i<no; i++) {
    xxx=100*(x1[i]-_min)/(_max-_min);
    yyy=100*(y1[i]-_ymin)/(_ymax-_ymin);
    if ((xxx>=0)&&(xxx<=100)&&(yyy>=0)&&(yyy<=100))hist.Put(round(xxx),round(yyy),0,hist.Get(round(xxx),round(yyy),0)+1);
  }
  _hist=hist;
  hist.Write("hist.nii.gz");

  //expectation values and covarinaces
  double x=0,y=0,xy=0,x2=0,y2=0;
  double ex,ey,ex2,exy,ey2;
  double varx, vary, covxy;
  double k,l,w;
  double sum=0;

  for (i=0; i<no; i++) {
    w=weights[i];
    sum +=w;

    x += x1[i]*w;
    y += y1[i]*w;
    xy += x1[i]*y1[i]*w;
    x2 += x1[i]*x1[i]*w;
    y2 += y1[i]*y1[i]*w;
  }

  ex  = x/sum;
  ey  = y/sum;
  ex2 = x2/sum;
  exy = xy/sum;
  ey2 = y2/sum;

  varx  = ex2-ex*ex;
  vary  = ey2-ey*ey;
  covxy = exy-ex*ey;

  k = (varx-vary)/(2*covxy);
  l = k+sqrt(1+k*k);

  _m = 1/l;
  _b = ey-ex/l;

  cerr<<"y = "<<_m<<" * x + "<<_b<<endl;

  double origx;
  for (i=0; i<=100; i++) {
    origx=_min+(_max-_min)*(i/100.0);
    yyy=100*(Lin(origx)-_ymin)/(_ymax-_ymin);
    if ((yyy>=0)&&(yyy<=100))trans.Put(i,round(yyy),0,1);
  }
  trans.Write("trans.nii.gz");
  _trans=trans;

}



