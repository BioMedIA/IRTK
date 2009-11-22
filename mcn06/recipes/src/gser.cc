#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7

void gser(float *gamser, float a, float x, float *gln)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	int n;
	float sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
