#include <math.h>

void simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp,
	float *bmax)
{
	int k;
	float test;

	*kp=ll[1];
	*bmax=a[mm+1][*kp+1];
	for (k=2;k<=nll;k++) {
		if (iabf == 0)
			test=a[mm+1][ll[k]+1]-(*bmax);
		else
			test=fabs(a[mm+1][ll[k]+1])-fabs(*bmax);
		if (test > 0.0) {
			*bmax=a[mm+1][ll[k]+1];
			*kp=ll[k];
		}
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
