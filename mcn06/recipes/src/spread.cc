#define NRANSI
#include "nrutil.h"

void spread(float y, float yy[], unsigned long n, float x, int m)
{
	int ihi,ilo,ix,j,nden;
	static int nfac[11]={0,1,1,2,6,24,120,720,5040,40320,362880};
	float fac;

	if (m > 10) nrerror("factorial table too small in spread");
	ix=(int)x;
	if (x == (float)ix) yy[ix] += y;
	else {
		ilo=LMIN(LMAX((long)(x-0.5*m+1.0),1),n-m+1);
		ihi=ilo+m-1;
		nden=nfac[m];
		fac=x-ilo;
		for (j=ilo+1;j<=ihi;j++) fac *= (x-j);
		yy[ihi] += y*fac/(nden*(x-ihi));
		for (j=ihi-1;j>=ilo;j--) {
			nden=(nden/(j+1-ilo))*(j-ihi);
			yy[j] += y*fac/(nden*(x-j));
		}
	}
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
