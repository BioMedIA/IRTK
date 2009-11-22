#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define NPFAC 8
#define MAXIT 5
#define PIO2 (3.141592653589793/2.0)
#define BIG 1.0e30

void ratlsq(double (*fn)(double), double a, double b, int mm, int kk,
	double cof[], double *dev)
{
	double ratval(double x, double cof[], int mm, int kk);
	void dsvbksb(double **u, double w[], double **v, int m, int n, double b[],
		double x[]);
	void dsvdcmp(double **a, int m, int n, double w[], double **v);
	int i,it,j,ncof,npt;
	double devmax,e,hth,power,sum,*bb,*coff,*ee,*fs,**u,**v,*w,*wt,*xs;

	ncof=mm+kk+1;
	npt=NPFAC*ncof;
	bb=dvector(1,npt);
	coff=dvector(0,ncof-1);
	ee=dvector(1,npt);
	fs=dvector(1,npt);
	u=dmatrix(1,npt,1,ncof);
	v=dmatrix(1,ncof,1,ncof);
	w=dvector(1,ncof);
	wt=dvector(1,npt);
	xs=dvector(1,npt);
	*dev=BIG;
	for (i=1;i<=npt;i++) {
		if (i < npt/2) {
			hth=PIO2*(i-1)/(npt-1.0);
			xs[i]=a+(b-a)*DSQR(sin(hth));
		} else {
			hth=PIO2*(npt-i)/(npt-1.0);
			xs[i]=b-(b-a)*DSQR(sin(hth));
		}
		fs[i]=(*fn)(xs[i]);
		wt[i]=1.0;
		ee[i]=1.0;
	}
	e=0.0;
	for (it=1;it<=MAXIT;it++) {
		for (i=1;i<=npt;i++) {
			power=wt[i];
			bb[i]=power*(fs[i]+SIGN(e,ee[i]));
			for (j=1;j<=mm+1;j++) {
				u[i][j]=power;
				power *= xs[i];
			}
			power = -bb[i];
			for (j=mm+2;j<=ncof;j++) {
				power *= xs[i];
				u[i][j]=power;
			}
		}
		dsvdcmp(u,npt,ncof,w,v);
		dsvbksb(u,w,v,npt,ncof,bb,coff-1);
		devmax=sum=0.0;
		for (j=1;j<=npt;j++) {
			ee[j]=ratval(xs[j],coff,mm,kk)-fs[j];
			wt[j]=fabs(ee[j]);
			sum += wt[j];
			if (wt[j] > devmax) devmax=wt[j];
		}
		e=sum/npt;
		if (devmax <= *dev) {
			for (j=0;j<ncof;j++) cof[j]=coff[j];
			*dev=devmax;
		}
		printf(" ratlsq iteration= %2d  max error= %10.3e\n",it,devmax);
	}
	free_dvector(xs,1,npt);
	free_dvector(wt,1,npt);
	free_dvector(w,1,ncof);
	free_dmatrix(v,1,ncof,1,ncof);
	free_dmatrix(u,1,npt,1,ncof);
	free_dvector(fs,1,npt);
	free_dvector(ee,1,npt);
	free_dvector(coff,0,ncof-1);
	free_dvector(bb,1,npt);
}
#undef NPFAC
#undef MAXIT
#undef PIO2
#undef BIG
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
