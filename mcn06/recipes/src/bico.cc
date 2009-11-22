#include <math.h>

float bico(int n, int k)
{
	float factln(int n);

	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
