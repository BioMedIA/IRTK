#define IAZ 48

void mp2dfr(unsigned char a[], unsigned char s[], int n, int *m)
{
	void mplsh(unsigned char u[], int n);
	void mpsmu(unsigned char w[], unsigned char u[], int n, int iv);
	int j;

	*m=(int) (2.408*n);
	for (j=1;j<=(*m);j++) {
		mpsmu(a,a,n,10);
		s[j]=a[1]+IAZ;
		mplsh(a,n);
	}
}
#undef IAZ
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
