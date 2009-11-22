#include <math.h>
#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-7

extern int ndatat;
extern float *xt,*yt,aa,abdevt;

float rofunc(float b)
{
	float select(unsigned long k, unsigned long n, float arr[]);
	int j;
	float *arr,d,sum=0.0;

	arr=vector(1,ndatat);
	for (j=1;j<=ndatat;j++) arr[j]=yt[j]-b*xt[j];
	if (ndatat & 1) {
		aa=select((ndatat+1)>>1,ndatat,arr);
	}
	else {
		j=ndatat >> 1;
		aa=0.5*(select(j,ndatat,arr)+select(j+1,ndatat,arr));
	}
	abdevt=0.0;
	for (j=1;j<=ndatat;j++) {
		d=yt[j]-(b*xt[j]+aa);
		abdevt += fabs(d);
		if (yt[j] != 0.0) d /= fabs(yt[j]);
		if (fabs(d) > EPS) sum += (d >= 0.0 ? xt[j] : -xt[j]);
	}
	free_vector(arr,1,ndatat);
	return sum;
}
#undef EPS
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
