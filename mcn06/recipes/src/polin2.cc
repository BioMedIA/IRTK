#define NRANSI
#include "nrutil.h"

void polin2(float x1a[], float x2a[], float **ya, int m, int n, float x1,
	float x2, float *y, float *dy)
{
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	int j;
	float *ymtmp;

	ymtmp=vector(1,m);
	for (j=1;j<=m;j++) {
		polint(x2a,ya[j],n,x2,&ymtmp[j],dy);
	}
	polint(x1a,ymtmp,m,x1,y,dy);
	free_vector(ymtmp,1,m);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
