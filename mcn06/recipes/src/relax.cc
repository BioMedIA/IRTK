void relax(double **u, double **rhs, int n)
{
	int i,ipass,isw,j,jsw=1;
	double h,h2;

	h=1.0/(n-1);
	h2=h*h;
	for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
		isw=jsw;
		for (j=2;j<n;j++,isw=3-isw)
			for (i=isw+1;i<n;i+=2)
				u[i][j]=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1]
					+u[i][j-1]-h2*rhs[i][j]);
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
