void copy(double **aout, double **ain, int n)
{
	int i,j;
	for (i=1;i<=n;i++)
		for (j=1;j<=n;j++)
			aout[j][i]=ain[j][i];

}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
