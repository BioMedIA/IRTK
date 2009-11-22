#include <math.h>
#define NRANSI
#include "nrutil.h"
#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

float dfridr(float (*func)(float), float x, float h, float *err)
{
	int i,j;
	float errt,fac,hh,**a,ans;

	if (h == 0.0) nrerror("h must be nonzero in dfridr.");
	a=matrix(1,NTAB,1,NTAB);
	hh=h;
	a[1][1]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
	*err=BIG;
	for (i=2;i<=NTAB;i++) {
		hh /= CON;
		a[1][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
		fac=CON2;
		for (j=2;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=FMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			if (errt <= *err) {
				*err=errt;
				ans=a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
	}
	free_matrix(a,1,NTAB,1,NTAB);
	return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
