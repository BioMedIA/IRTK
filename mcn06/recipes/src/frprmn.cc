#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []))
{
	void linmin(float p[], float xi[], int n, float *fret,
		float (*func)(float []));
	int j,its;
	float gg,gam,fp,dgg;
	float *g,*h,*xi;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		linmin(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
