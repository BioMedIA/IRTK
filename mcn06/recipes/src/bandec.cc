#include <math.h>
#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20

void bandec(float **a, unsigned long n, int m1, int m2, float **al,
	unsigned long indx[], float *d)
{
	unsigned long i,j,k,l;
	int mm;
	float dum;

	mm=m1+m2+1;
	l=m1;
	for (i=1;i<=m1;i++) {
		for (j=m1+2-i;j<=mm;j++) a[i][j-l]=a[i][j];
		l--;
		for (j=mm-l;j<=mm;j++) a[i][j]=0.0;
	}
	*d=1.0;
	l=m1;
	for (k=1;k<=n;k++) {
		dum=a[k][1];
		i=k;
		if (l < n) l++;
		for (j=k+1;j<=l;j++) {
			if (fabs(a[j][1]) > fabs(dum)) {
				dum=a[j][1];
				i=j;
			}
		}
		indx[k]=i;
		if (dum == 0.0) a[k][1]=TINY;
		if (i != k) {
			*d = -(*d);
			for (j=1;j<=mm;j++) SWAP(a[k][j],a[i][j])
		}
		for (i=k+1;i<=l;i++) {
			dum=a[i][1]/a[k][1];
			al[k][i-k]=dum;
			for (j=2;j<=mm;j++) a[i][j-1]=a[i][j]-dum*a[k][j];
			a[i][mm]=0.0;
		}
	}
}
#undef SWAP
#undef TINY
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
