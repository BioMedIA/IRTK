#include <math.h>
#define MAXITS 1000
#define EPS 1.0e-5

void sor(double **a, double **b, double **c, double **d, double **e,
	double **f, double **u, int jmax, double rjac)
{
	void nrerror(char error_text[]);
	int ipass,j,jsw,l,lsw,n;
	double anorm,anormf=0.0,omega=1.0,resid;

	for (j=2;j<jmax;j++)
		for (l=2;l<jmax;l++)
			anormf += fabs(f[j][l]);
	for (n=1;n<=MAXITS;n++) {
		anorm=0.0;
		jsw=1;
		for (ipass=1;ipass<=2;ipass++) {
			lsw=jsw;
			for (j=2;j<jmax;j++) {
				for (l=lsw+1;l<jmax;l+=2) {
					resid=a[j][l]*u[j+1][l]
						+b[j][l]*u[j-1][l]
						+c[j][l]*u[j][l+1]
						+d[j][l]*u[j][l-1]
						+e[j][l]*u[j][l]
						-f[j][l];
					anorm += fabs(resid);
					u[j][l] -= omega*resid/e[j][l];
				}
				lsw=3-lsw;
			}
			jsw=3-jsw;
			omega=(n == 1 && ipass == 1 ? 1.0/(1.0-0.5*rjac*rjac) :
				1.0/(1.0-0.25*rjac*rjac*omega));
		}
		if (anorm < EPS*anormf) return;
	}
	nrerror("MAXITS exceeded");
}
#undef MAXITS
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
