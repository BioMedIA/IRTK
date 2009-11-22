void addint(double **uf, double **uc, double **res, int nf)
{
	void interp(double **uf, double **uc, int nf);
	int i,j;

	interp(res,uc,nf);
	for (j=1;j<=nf;j++)
		for (i=1;i<=nf;i++)
			uf[i][j] += res[i][j];
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
