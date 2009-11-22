#include <math.h>
#define NRANSI
#include "nrutil.h"

void rotate(float **r, float **qt, int n, int i, float a, float b)
{
	int j;
	float c,fact,s,w,y;

	if (a == 0.0) {
		c=0.0;
		s=(b >= 0.0 ? 1.0 : -1.0);
	} else if (fabs(a) > fabs(b)) {
		fact=b/a;
		c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
		s=fact*c;
	} else {
		fact=a/b;
		s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
		c=fact*s;
	}
	for (j=i;j<=n;j++) {
		y=r[i][j];
		w=r[i+1][j];
		r[i][j]=c*y-s*w;
		r[i+1][j]=s*y+c*w;
	}
	for (j=1;j<=n;j++) {
		y=qt[i][j];
		w=qt[i+1][j];
		qt[i][j]=c*y-s*w;
		qt[i+1][j]=s*y+c*w;
	}
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
