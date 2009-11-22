#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}

void banbks(float **a, unsigned long n, int m1, int m2, float **al,
	unsigned long indx[], float b[])
{
	unsigned long i,k,l;
	int mm;
	float dum;

	mm=m1+m2+1;
	l=m1;
	for (k=1;k<=n;k++) {
		i=indx[k];
		if (i != k) SWAP(b[k],b[i])
		if (l < n) l++;
		for (i=k+1;i<=l;i++) b[i] -= al[k][i-k]*b[k];
	}
	l=1;
	for (i=n;i>=1;i--) {
		dum=b[i];
		for (k=2;k<=l;k++) dum -= a[i][k]*b[k+i-1];
		b[i]=dum/a[i][1];
		if (l < mm) l++;
	}
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
