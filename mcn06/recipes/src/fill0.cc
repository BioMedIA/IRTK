void fill0(double **u, int n)
{
	int i,j;
	for (j=1;j<=n;j++)
		for (i=1;i<=n;i++)
			u[i][j]=0.0;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
