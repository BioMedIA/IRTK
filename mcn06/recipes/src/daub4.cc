#define NRANSI
#include "nrutil.h"
#define C0 0.4829629131445341
#define C1 0.8365163037378079
#define C2 0.2241438680420134
#define C3 -0.1294095225512604

void daub4(float a[], unsigned long n, int isign)
{
	float *wksp;
	unsigned long nh,nh1,i,j;

	if (n < 4) return;
	wksp=vector(1,n);
	nh1=(nh=n >> 1)+1;
	if (isign >= 0) {
		for (i=1,j=1;j<=n-3;j+=2,i++) {
			wksp[i]=C0*a[j]+C1*a[j+1]+C2*a[j+2]+C3*a[j+3];
			wksp[i+nh] = C3*a[j]-C2*a[j+1]+C1*a[j+2]-C0*a[j+3];
		}
		wksp[i]=C0*a[n-1]+C1*a[n]+C2*a[1]+C3*a[2];
		wksp[i+nh] = C3*a[n-1]-C2*a[n]+C1*a[1]-C0*a[2];
	} else {
		wksp[1]=C2*a[nh]+C1*a[n]+C0*a[1]+C3*a[nh1];
		wksp[2] = C3*a[nh]-C0*a[n]+C1*a[1]-C2*a[nh1];
		for (i=1,j=3;i<nh;i++) {
			wksp[j++]=C2*a[i]+C1*a[i+nh]+C0*a[i+1]+C3*a[i+nh1];
			wksp[j++] = C3*a[i]-C0*a[i+nh]+C1*a[i+1]-C2*a[i+nh1];
		}
	}
	for (i=1;i<=n;i++) a[i]=wksp[i];
	free_vector(wksp,1,n);
}
#undef C0
#undef C1
#undef C2
#undef C3
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
