void matsub(double **a, double **b, double **c, int n)
{
	int i,j;

	for (j=1;j<=n;j++)
		for (i=1;i<=n;i++)
			c[i][j]=a[i][j]-b[i][j];
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
