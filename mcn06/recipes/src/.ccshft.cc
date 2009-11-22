void pcshft(float a, float b, float d[], int n)
{
	int k,j;
	float fac,cnst;

	cnst=2.0/(b-a);
	fac=cnst;
	for (j=1;j<n;j++) {
		d[j] *= fac;
		fac *= cnst;
	}
	cnst=0.5*(a+b);
	for (j=0;j<=n-2;j++)
		for (k=n-2;k>=j;k--)
			d[k] -= cnst*d[k+1];
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
