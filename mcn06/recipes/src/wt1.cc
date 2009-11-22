void wt1(float a[], unsigned long n, int isign,
	void (*wtstep)(float [], unsigned long, int))
{
	unsigned long nn;

	if (n < 4) return;
	if (isign >= 0) {
		for (nn=n;nn>=4;nn>>=1) (*wtstep)(a,nn,isign);
	} else {
		for (nn=4;nn<=n;nn<<=1) (*wtstep)(a,nn,isign);
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
