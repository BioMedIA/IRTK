#include <math.h>

void fgauss(float x, float a[], float *y, float dyda[], int na)
{
	int i;
	float fac,ex,arg;

	*y=0.0;
	for (i=1;i<=na-1;i+=3) {
		arg=(x-a[i+1])/a[i+2];
		ex=exp(-arg*arg);
		fac=a[i]*ex*2.0*arg;
		*y += a[i]*ex;
		dyda[i]=ex;
		dyda[i+1]=fac/a[i+2];
		dyda[i+2]=fac*arg/a[i+2];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
