#include <math.h>
#define NRANSI
#include "nrutil.h"

float ellpi(float phi, float en, float ak)
{
	float rf(float x, float y, float z);
	float rj(float x, float y, float z, float p);
	float cc,enss,q,s;

	s=sin(phi);
	enss=en*s*s;
	cc=SQR(cos(phi));
	q=(1.0-s*ak)*(1.0+s*ak);
	return s*(rf(cc,q,1.0)-enss*rj(cc,q,1.0,1.0+enss)/3.0);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
