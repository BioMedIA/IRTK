#define NRANSI
#include "nrutil.h"

typedef struct {
	int ncof,ioff,joff;
	float *cc,*cr;
} wavefilt;

extern wavefilt wfilt;

void pwt(float a[], unsigned long n, int isign)
{
	float ai,ai1,*wksp;
	unsigned long i,ii,j,jf,jr,k,n1,ni,nj,nh,nmod;

	if (n < 4) return;
	wksp=vector(1,n);
	nmod=wfilt.ncof*n;
	n1=n-1;
	nh=n >> 1;
	for (j=1;j<=n;j++) wksp[j]=0.0;
	if (isign >= 0) {
		for (ii=1,i=1;i<=n;i+=2,ii++) {
			ni=i+nmod+wfilt.ioff;
			nj=i+nmod+wfilt.joff;
			for (k=1;k<=wfilt.ncof;k++) {
				jf=n1 & (ni+k);
				jr=n1 & (nj+k);
				wksp[ii] += wfilt.cc[k]*a[jf+1];
				wksp[ii+nh] += wfilt.cr[k]*a[jr+1];
			}
		}
	} else {
		for (ii=1,i=1;i<=n;i+=2,ii++) {
			ai=a[ii];
			ai1=a[ii+nh];
			ni=i+nmod+wfilt.ioff;
			nj=i+nmod+wfilt.joff;
			for (k=1;k<=wfilt.ncof;k++) {
				jf=(n1 & (ni+k))+1;
				jr=(n1 & (nj+k))+1;
				wksp[jf] += wfilt.cc[k]*ai;
				wksp[jr] += wfilt.cr[k]*ai1;
			}
		}
	}
	for (j=1;j<=n;j++) a[j]=wksp[j];
	free_vector(wksp,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
