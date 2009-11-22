#include <math.h>
#define NRANSI
#include "nrutil.h"
#define TWOPID 6.2831853071795865

void period(float x[], float y[], int n, float ofac, float hifac, float px[],
	float py[], int np, int *nout, int *jmax, float *prob)
{
	void avevar(float data[], unsigned long n, float *ave, float *var);
	int i,j;
	float ave,c,cc,cwtau,effm,expy,pnow,pymax,s,ss,sumc,sumcy,sums,sumsh,
		sumsy,swtau,var,wtau,xave,xdif,xmax,xmin,yy;
	double arg,wtemp,*wi,*wpi,*wpr,*wr;

	wi=dvector(1,n);
	wpi=dvector(1,n);
	wpr=dvector(1,n);
	wr=dvector(1,n);
	*nout=0.5*ofac*hifac*n;
	if (*nout > np) nrerror("output arrays too short in period");
	avevar(y,n,&ave,&var);
	xmax=xmin=x[1];
	for (j=1;j<=n;j++) {
		if (x[j] > xmax) xmax=x[j];
		if (x[j] < xmin) xmin=x[j];
	}
	xdif=xmax-xmin;
	xave=0.5*(xmax+xmin);
	pymax=0.0;
	pnow=1.0/(xdif*ofac);
	for (j=1;j<=n;j++) {
		arg=TWOPID*((x[j]-xave)*pnow);
		wpr[j] = -2.0*SQR(sin(0.5*arg));
		wpi[j]=sin(arg);
		wr[j]=cos(arg);
		wi[j]=wpi[j];
	}
	for (i=1;i<=(*nout);i++) {
		px[i]=pnow;
		sumsh=sumc=0.0;
		for (j=1;j<=n;j++) {
			c=wr[j];
			s=wi[j];
			sumsh += s*c;
			sumc += (c-s)*(c+s);
		}
#ifdef WIN32
    wtau=0.5*atan2(static_cast<double>(2.0*sumsh),static_cast<double>(sumc));
#else
		wtau=0.5*atan2(2.0*sumsh,sumc);
#endif
		swtau=sin(wtau);
		cwtau=cos(wtau);
		sums=sumc=sumsy=sumcy=0.0;
		for (j=1;j<=n;j++) {
			s=wi[j];
			c=wr[j];
			ss=s*cwtau-c*swtau;
			cc=c*cwtau+s*swtau;
			sums += ss*ss;
			sumc += cc*cc;
			yy=y[j]-ave;
			sumsy += yy*ss;
			sumcy += yy*cc;
			wr[j]=((wtemp=wr[j])*wpr[j]-wi[j]*wpi[j])+wr[j];
			wi[j]=(wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
		}
		py[i]=0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;
		if (py[i] >= pymax) pymax=py[(*jmax=i)];
		pnow += 1.0/(ofac*xdif);
	}
	expy=exp(-pymax);
	effm=2.0*(*nout)/ofac;
	*prob=effm*expy;
	if (*prob > 0.01)
#ifdef WIN32
    *prob=1.0-pow(static_cast<double>(1.0-expy),static_cast<double>(effm));
#else
    *prob=1.0-pow(1.0-expy,effm);
#endif
	free_dvector(wr,1,n);
	free_dvector(wpr,1,n);
	free_dvector(wpi,1,n);
	free_dvector(wi,1,n);
}
#undef TWOPID
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
