#include <math.h>

float expdev(long *idum)
{
	float ran1(long *idum);
	float dum;

	do
		dum=ran1(idum);
	while (dum == 0.0);
	return -log(dum);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
