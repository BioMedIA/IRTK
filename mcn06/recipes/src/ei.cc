#include <math.h>
#define EULER 0.57721566
#define MAXIT 100
#define FPMIN 1.0e-30
#define EPS 6.0e-8

float ei(float x)
{
	void nrerror(char error_text[]);
	int k;
	float fact,prev,sum,term;

	if (x <= 0.0) nrerror("Bad argument in ei");
	if (x < FPMIN) return log(x)+EULER;
	if (x <= -log(EPS)) {
		sum=0.0;
		fact=1.0;
		for (k=1;k<=MAXIT;k++) {
			fact *= x/k;
			term=fact/k;
			sum += term;
			if (term < EPS*sum) break;
		}
		if (k > MAXIT) nrerror("Series failed in ei");
		return sum+log(x)+EULER;
	} else {
		sum=0.0;
		term=1.0;
		for (k=1;k<=MAXIT;k++) {
			prev=term;
			term *= k/x;
			if (term < EPS) break;
			if (term < prev) sum += term;
			else {
				sum -= prev;
				break;
			}
		}
		return exp(x)*(1.0+sum)/x;
	}
}
#undef EPS
#undef EULER
#undef MAXIT
#undef FPMIN
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
