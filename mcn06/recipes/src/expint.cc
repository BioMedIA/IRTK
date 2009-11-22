#include <math.h>
#define MAXIT 100
#define EULER 0.5772156649
#define FPMIN 1.0e-30
#define EPS 1.0e-7

float expint(int n, float x)
{
	void nrerror(char error_text[]);
	int i,ii,nm1;
	float a,b,c,d,del,fact,h,psi,ans;

	nm1=n-1;
	if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1)))
	nrerror("bad arguments in expint");
	else {
		if (n == 0) ans=exp(-x)/x;
		else {
			if (x == 0.0) ans=1.0/nm1;

			else {
				if (x > 1.0) {
					b=x+n;
					c=1.0/FPMIN;
					d=1.0/b;
					h=d;
					for (i=1;i<=MAXIT;i++) {
						a = -i*(nm1+i);
						b += 2.0;
						d=1.0/(a*d+b);
						c=b+a/c;
						del=c*d;
						h *= del;
						if (fabs(del-1.0) < EPS) {
							ans=h*exp(-x);
							return ans;
						}
					}
					nrerror("continued fraction failed in expint");
				} else {
					ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);
					fact=1.0;
					for (i=1;i<=MAXIT;i++) {
						fact *= -x/i;
						if (i != nm1) del = -fact/(i-nm1);
						else {
							psi = -EULER;
							for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
							del=fact*(-log(x)+psi);
						}
						ans += del;
						if (fabs(del) < fabs(ans)*EPS) return ans;
					}
					nrerror("series failed in expint");
				}
			}
		}
	}
	return ans;
}
#undef MAXIT
#undef EPS
#undef FPMIN
#undef EULER
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
