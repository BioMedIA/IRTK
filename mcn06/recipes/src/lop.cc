void lop(double **out, double **u, int n)
{
	int i,j;
	double h,h2i;

	h=1.0/(n-1);
	h2i=1.0/(h*h);
	for (j=2;j<n;j++)
		for (i=2;i<n;i++)
			out[i][j]=h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
				4.0*u[i][j])+u[i][j]*u[i][j];
	for (i=1;i<=n;i++)
		out[i][1]=out[i][n]=out[1][i]=out[n][i]=0.0;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
