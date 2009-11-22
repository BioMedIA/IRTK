#include <math.h>
#define RAD (3.14159265/180.0)

void flmoon(int n, int nph, long *jd, float *frac)
{
	void nrerror(char error_text[]);
	int i;
	float am,as,c,t,t2,xtra;

	c=n+nph/4.0;
	t=c/1236.85;
	t2=t*t;
	as=359.2242+29.105356*c;
	am=306.0253+385.816918*c+0.010730*t2;
	*jd=2415020+28L*n+7L*nph;
	xtra=0.75933+1.53058868*c+((1.178e-4)-(1.55e-7)*t)*t2;
	if (nph == 0 || nph == 2)
		xtra += (0.1734-3.93e-4*t)*sin(RAD*as)-0.4068*sin(RAD*am);
	else if (nph == 1 || nph == 3)
		xtra += (0.1721-4.0e-4*t)*sin(RAD*as)-0.6280*sin(RAD*am);
	else nrerror("nph is unknown in flmoon");
	i=(int)(xtra >= 0.0 ? floor(xtra) : ceil(xtra-1.0));
	*jd += i;
	*frac=xtra-i;
}
#undef RAD
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
