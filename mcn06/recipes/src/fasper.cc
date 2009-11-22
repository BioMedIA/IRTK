#include <math.h>
#define NRANSI
#include "nrutil.h"
#define MOD(a,b) 	while(a >= b) a -= b;
#define MACC 4

void fasper(float x[], float y[], unsigned long n, float ofac, float hifac,
	float wk1[], float wk2[], unsigned long nwk, unsigned long *nout,
	unsigned long *jmax, float *prob)
{
	void avevar(float data[], unsigned long n, float *ave, float *var);
	void realft(float data[], unsigned long n, int isign);
	void spread(float y, float yy[], unsigned long n, float x, int m);
	unsigned long j,k,ndim,nfreq,nfreqt;
	float ave,ck,ckk,cterm,cwt,den,df,effm,expy,fac,fndim,hc2wt;
	float hs2wt,hypo,pmax,sterm,swt,var,xdif,xmax,xmin;

	*nout=0.5*ofac*hifac*n;
	nfreqt=ofac*hifac*n*MACC;
	nfreq=64;
	while (nfreq < nfreqt) nfreq <<= 1;
	ndim=nfreq << 1;
	if (ndim > nwk) nrerror("workspaces too small in fasper");
	avevar(y,n,&ave,&var);
	xmin=x[1];
	xmax=xmin;
	for (j=2;j<=n;j++) {
		if (x[j] < xmin) xmin=x[j];
		if (x[j] > xmax) xmax=x[j];
	}
	xdif=xmax-xmin;
	for (j=1;j<=ndim;j++) wk1[j]=wk2[j]=0.0;
	fac=ndim/(xdif*ofac);
	fndim=ndim;
	for (j=1;j<=n;j++) {
		ck=(x[j]-xmin)*fac;
		MOD(ck,fndim)
		ckk=2.0*(ck++);
		MOD(ckk,fndim)
		++ckk;
		spread(y[j]-ave,wk1,ndim,ck,MACC);
		spread(1.0,wk2,ndim,ckk,MACC);
	}
	realft(wk1,ndim,1);
	realft(wk2,ndim,1);
	df=1.0/(xdif*ofac);
	pmax = -1.0;
	for (k=3,j=1;j<=(*nout);j++,k+=2) {
		hypo=sqrt(wk2[k]*wk2[k]+wk2[k+1]*wk2[k+1]);
		hc2wt=0.5*wk2[k]/hypo;
		hs2wt=0.5*wk2[k+1]/hypo;
		cwt=sqrt(0.5+hc2wt);
		swt=SIGN(sqrt(0.5-hc2wt),hs2wt);
		den=0.5*n+hc2wt*wk2[k]+hs2wt*wk2[k+1];
		cterm=SQR(cwt*wk1[k]+swt*wk1[k+1])/den;
		sterm=SQR(cwt*wk1[k+1]-swt*wk1[k])/(n-den);
		wk1[j]=j*df;
		wk2[j]=(cterm+sterm)/(2.0*var);
		if (wk2[j] > pmax) pmax=wk2[(*jmax=j)];
	}
	expy=exp(-pmax);
	effm=2.0*(*nout)/ofac;
	*prob=effm*expy;
	if (*prob > 0.01)
#ifdef WIN32
    *prob=1.0-pow(static_cast<double>(1.0-expy),static_cast<double>(effm));
#else
    *prob=1.0-pow(1.0-expy,effm);
#endif
}
#undef MOD
#undef MACC
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
