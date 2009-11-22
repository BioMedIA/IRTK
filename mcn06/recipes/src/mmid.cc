#define NRANSI
#include "nrutil.h"

void mmid(float y[], float dydx[], int nvar, float xs, float htot, int nstep,
	float yout[], void (*derivs)(float, float[], float[]))
{
	int n,i;
	float x,swap,h2,h,*ym,*yn;

	ym=vector(1,nvar);
	yn=vector(1,nvar);
	h=htot/nstep;
	for (i=1;i<=nvar;i++) {
		ym[i]=y[i];
		yn[i]=y[i]+h*dydx[i];
	}
	x=xs+h;
	(*derivs)(x,yn,yout);
	h2=2.0*h;
	for (n=2;n<=nstep;n++) {
		for (i=1;i<=nvar;i++) {
			swap=ym[i]+h2*yout[i];
			ym[i]=yn[i];
			yn[i]=swap;
		}
		x += h;
		(*derivs)(x,yn,yout);
	}
	for (i=1;i<=nvar;i++)
		yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
	free_vector(yn,1,nvar);
	free_vector(ym,1,nvar);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
