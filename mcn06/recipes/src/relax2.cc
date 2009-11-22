void relax2(double **u, double **rhs, int n)
{
	int i,ipass,isw,j,jsw=1;
	double foh2,h,h2i,res;

	h=1.0/(n-1);
	h2i=1.0/(h*h);
	foh2 = -4.0*h2i;
	for (ipass=1;ipass<=2;ipass++,jsw=3-jsw) {
		isw=jsw;
		for (j=2;j<n;j++,isw=3-isw) {
			for (i=isw+1;i<n;i+=2) {
				res=h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
					4.0*u[i][j])+u[i][j]*u[i][j]-rhs[i][j];
				u[i][j] -= res/(foh2+2.0*u[i][j]);
			}
		}
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
