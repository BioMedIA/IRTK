#define NRANSI
#include "nrutil.h"

void vander(double x[], double w[], double q[], int n)
{
	int i,j,k,k1;
	double b,s,t,xx;
	double *c;

	c=dvector(1,n);
	if (n == 1) w[1]=q[1];
	else {
		for (i=1;i<=n;i++) c[i]=0.0;
		c[n] = -x[1];
		for (i=2;i<=n;i++) {
			xx = -x[i];
			for (j=(n+1-i);j<=(n-1);j++) c[j] += xx*c[j+1];
			c[n] += xx;
		}
		for (i=1;i<=n;i++) {
			xx=x[i];
			t=b=1.0;
			s=q[n];
			k=n;
			for (j=2;j<=n;j++) {
				k1=k-1;
				b=c[k]+xx*b;
				s += q[k1]*b;
				t=xx*t+b;
				k=k1;
			}
			w[i]=s/t;
		}
	}
	free_dvector(c,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
