#include <math.h>

float beta(float z, float w)
{
	float gammln(float xx);

	return exp(gammln(z)+gammln(w)-gammln(z+w));
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
