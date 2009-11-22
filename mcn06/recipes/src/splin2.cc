#define NRANSI
#include "nrutil.h"

void splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n,
	float x1, float x2, float *y)
{
	void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
	void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
	int j;
	float *ytmp,*yytmp;

	ytmp=vector(1,m);
	yytmp=vector(1,m);
	for (j=1;j<=m;j++)
		splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
	spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
	splint(x1a,yytmp,ytmp,m,x1,y);
	free_vector(yytmp,1,m);
	free_vector(ytmp,1,m);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
