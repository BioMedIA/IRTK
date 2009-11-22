#include <math.h>
#define NRANSI
#include "nrutil.h"

float elle(float phi, float ak)
{
	float rd(float x, float y, float z);
	float rf(float x, float y, float z);
	float cc,q,s;

	s=sin(phi);
	cc=SQR(cos(phi));
	q=(1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-(SQR(s*ak))*rd(cc,q,1.0)/3.0);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
