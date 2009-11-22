#include <math.h>
#define NRANSI
#include "nrutil.h"

void qrdcmp(float **a, int n, float *c, float *d, int *sing)
{
	int i,j,k;
	float scale=0.0,sigma,sum,tau;

	*sing=0;
	for (k=1;k<n;k++) {
		for (i=k;i<=n;i++) scale=FMAX(scale,fabs(a[i][k]));
		if (scale == 0.0) {
			*sing=1;
			c[k]=d[k]=0.0;
		} else {
			for (i=k;i<=n;i++) a[i][k] /= scale;
			for (sum=0.0,i=k;i<=n;i++) sum += SQR(a[i][k]);
			sigma=SIGN(sqrt(sum),a[k][k]);
			a[k][k] += sigma;
			c[k]=sigma*a[k][k];
			d[k] = -scale*sigma;
			for (j=k+1;j<=n;j++) {
				for (sum=0.0,i=k;i<=n;i++) sum += a[i][k]*a[i][j];
				tau=sum/c[k];
				for (i=k;i<=n;i++) a[i][j] -= tau*a[i][k];
			}
		}
	}
	d[n]=a[n][n];
	if (d[n] == 0.0) *sing=1;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
