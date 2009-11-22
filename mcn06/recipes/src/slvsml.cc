void slvsml(double **u, double **rhs)
{
	void fill0(double **u, int n);
	double h=0.5;

	fill0(u,3);
	u[2][2] = -h*h*rhs[2][2]/4.0;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
