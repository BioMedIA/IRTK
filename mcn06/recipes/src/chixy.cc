#include <math.h>
#define NRANSI
#include "nrutil.h"
#define BIG 1.0e30

extern int nn;
extern float *xx,*yy,*sx,*sy,*ww,aa,offs;

float chixy(float bang)
{
	int j;
	float ans,avex=0.0,avey=0.0,sumw=0.0,b;

	b=tan(bang);
	for (j=1;j<=nn;j++) {
		ww[j] = SQR(b*sx[j])+SQR(sy[j]);
		sumw += (ww[j] = (ww[j] == 0.0 ? BIG : 1.0/ww[j]));
		avex += ww[j]*xx[j];
		avey += ww[j]*yy[j];
	}
	if (sumw == 0.0) sumw = BIG;
	avex /= sumw;
	avey /= sumw;
	aa=avey-b*avex;
	for (ans = -offs,j=1;j<=nn;j++)
		ans += ww[j]*SQR(yy[j]-aa-b*xx[j]);
	return ans;
}
#undef BIG
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
