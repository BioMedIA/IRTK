#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n); \
free_matrix(hessin,1,n,1,n);free_vector(hdg,1,n);free_vector(g,1,n); \
free_vector(dg,1,n);

void dfpmin(float p[], int n, float gtol, int *iter, float *fret,
	float(*func)(float []), void (*dfunc)(float [], float []))
{
	void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
		 float *f, float stpmax, int *check, float (*func)(float []));
	int check,i,its,j;
	float den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
	float *dg,*g,*hdg,**hessin,*pnew,*xi;

	dg=vector(1,n);
	g=vector(1,n);
	hdg=vector(1,n);
	hessin=matrix(1,n,1,n);
	pnew=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,g);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*FMAX(sqrt(sum),(float)n);
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
		fp = *fret;
		for (i=1;i<=n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i];
		(*dfunc)(p,g);
		test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
			FREEALL
			return;
		}
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac*fac > EPS*sumdg*sumxi) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
				}
			}
		}
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	nrerror("too many iterations in dfpmin");
	FREEALL
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
