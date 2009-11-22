#include <math.h>

float factrl(int n)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	static int ntop=4;
	static float a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;

	if (n < 0) nrerror("Negative factorial in routine factrl");
	if (n > 32) return exp(gammln(n+1.0));
	while (ntop<n) {
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
