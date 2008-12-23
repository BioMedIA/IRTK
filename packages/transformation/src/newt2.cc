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

#include <irtkTransformation.h>

#define NRANSI
#include <nrutil.h>
#define MAXITS 100
#define TOLF 1.0e-4
#define ALF 1.0e-4
#define TOLMIN 1.0e-6
#define TOLX 1.0e-7
#define STPMX 100.0

irtkTransformation *irtkTransformationPointer;

int nn;
float *fvec;
void (*nrfuncv)(int, float v[], float f[]);
double x_invert, y_invert, z_invert;

void newt2(float x[], int n, int *check, void (*vecfunc)(int, float [], float []))
{
  void fdjac(int n, float x[], float fvec[], float **df,
             void (*vecfunc)(int, float [], float []));
  float fmin(float x[]);
  void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
              float *f, float stpmax, int *check, float (*func)(float []));
  void lubksb(float **a, int n, int *indx, float b[]);
  void ludcmp(float **a, int n, int *indx, float *d);
  int i,its,j,*indx;
  float d,den,f,fold,stpmax,sum,temp,test,**fjac,*g,*p,*xold;
  double xpoint,ypoint,zpoint;
  irtkMatrix jacobian(3,3); /* Jacobian of transformation*/
  indx=ivector(1,n);
  fjac=matrix(1,n,1,n);
  g=vector(1,n);
  p=vector(1,n);
  xold=vector(1,n);

  fvec=vector(1,n);
  nn=n;
  nrfuncv=vecfunc;
  f=fmin(x);
  test=0.0;
  for (i=1;i<=n;i++)
    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
  if (test<0.01*TOLF) {
    *check=0;
    free_vector(fvec,1,n);
    free_vector(xold,1,n);
    free_vector(p,1,n);
    free_vector(g,1,n);
    free_matrix(fjac,1,n,1,n);
    free_ivector(indx,1,n);
    return;
  }
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
  stpmax=STPMX*FMAX(sqrt(sum),(float)n);
  for (its=1;its<=MAXITS;its++) {
    xpoint=x[1];
    ypoint=x[2];
    zpoint=x[3];
    irtkTransformationPointer->Jacobian(jacobian,xpoint,ypoint,zpoint);
    if (irtkTransformationPointer->Jacobian(xpoint,ypoint,zpoint) < 0.0) {
      cerr << "input mapping non-invertible!\n";
    }
    for (i=1;i<=n;i++) {
      for (j=1;j<=n;j++) fjac[j][i]=jacobian(j-1,i-1);
    }

    for (i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
      g[i]=sum;
    }
    for (i=1;i<=n;i++) xold[i]=x[i];
    fold=f;
    for (i=1;i<=n;i++) p[i] = -fvec[i];
    ludcmp(fjac,n,indx,&d);

    lubksb(fjac,n,indx,p);

    lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
    test=0.0;
    for (i=1;i<=n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) {
      *check=0;
      free_vector(fvec,1,n);
      free_vector(xold,1,n);
      free_vector(p,1,n);
      free_vector(g,1,n);
      free_matrix(fjac,1,n,1,n);
      free_ivector(indx,1,n);
      return;
    }

    if (*check) {
      test=0.0;
      den=FMAX(f,0.5*n);
      for (i=1;i<=n;i++) {
        temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
        if (temp > test) test=temp;
      }
      *check=(test < TOLMIN ? 1 : 0);
      free_vector(fvec,1,n);
      free_vector(xold,1,n);
      free_vector(p,1,n);
      free_vector(g,1,n);
      free_matrix(fjac,1,n,1,n);
      free_ivector(indx,1,n);
      return;
    }
    test=0.0;
    for (i=1;i<=n;i++) {
      temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) {
      free_vector(fvec,1,n);
      free_vector(xold,1,n);
      free_vector(p,1,n);
      free_vector(g,1,n);
      free_matrix(fjac,1,n,1,n);
      free_ivector(indx,1,n);
      return;
    }
  }
  *check=2;
  free_vector(fvec,1,n);
  free_vector(xold,1,n);
  free_vector(p,1,n);
  free_vector(g,1,n);
  free_matrix(fjac,1,n,1,n);
  free_ivector(indx,1,n);
  return;
}

void irtkTransformationEvaluate(int, float point[], float fval[])
{
  double xpoint,ypoint,zpoint;

  xpoint = point[1];
  ypoint = point[2];
  zpoint = point[3];

  irtkTransformationPointer->Transform(xpoint, ypoint, zpoint);

  fval[1] = xpoint-x_invert;
  fval[2] = ypoint-y_invert;
  fval[3] = zpoint-z_invert;
}

