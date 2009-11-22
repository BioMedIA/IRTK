#include <math.h>
#define JMAX 40

float rtbis(float (*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float dx,f,fmid,xmid,rtb;

	f=(*func)(x1);
	fmid=(*func)(x2);
	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}
#undef JMAX
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
