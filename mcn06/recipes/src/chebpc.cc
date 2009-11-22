#define NRANSI
#include "nrutil.h"

void chebpc(float c[], float d[], int n)
{
	int k,j;
	float sv,*dd;

	dd=vector(0,n-1);
	for (j=0;j<n;j++) d[j]=dd[j]=0.0;
	d[0]=c[n-1];
	for (j=n-2;j>=1;j--) {
		for (k=n-j;k>=1;k--) {
			sv=d[k];
			d[k]=2.0*d[k-1]-dd[k];
			dd[k]=sv;
		}
		sv=d[0];
		d[0] = -dd[0]+c[j];
		dd[0]=sv;
	}
	for (j=n-1;j>=1;j--)
		d[j]=d[j-1]-dd[j];
	d[0] = -dd[0]+0.5*c[0];
	free_vector(dd,0,n-1);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
