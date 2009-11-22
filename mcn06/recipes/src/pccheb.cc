void pccheb(float d[], float c[], int n)
{
	int j,jm,jp,k;
	float fac,pow;

	pow=1.0;
	c[0]=2.0*d[0];
	for (k=1;k<n;k++) {
		c[k]=0.0;
		fac=d[k]/pow;
		jm=k;
		jp=1;
		for (j=k;j>=0;j-=2,jm--,jp++) {
			c[j] += fac;
			fac *= ((float)jm)/((float)jp);
		}
		pow += pow;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
