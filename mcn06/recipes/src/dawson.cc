#include <math.h>
#define NRANSI
#include "nrutil.h"
#define NMAX 6
#define H 0.4
#define A1 (2.0/3.0)
#define A2 0.4
#define A3 (2.0/7.0)

float dawson(float x)
{
	int i,n0;
	float d1,d2,e1,e2,sum,x2,xp,xx,ans;
	static float c[NMAX+1];
	static int init = 0;

	if (init == 0) {
		init=1;
		for (i=1;i<=NMAX;i++) c[i]=exp(-SQR((2.0*i-1.0)*H));
	}
	if (fabs(x) < 0.2) {
		x2=x*x;
		ans=x*(1.0-A1*x2*(1.0-A2*x2*(1.0-A3*x2)));
	} else {
		xx=fabs(x);
		n0=2*(int)(0.5*xx/H+0.5);
		xp=xx-n0*H;
		e1=exp(2.0*xp*H);
		e2=e1*e1;
		d1=n0+1;
		d2=d1-2.0;
		sum=0.0;
		for (i=1;i<=NMAX;i++,d1+=2.0,d2-=2.0,e1*=e2)
			sum += c[i]*(e1/d1+1.0/(d2*e1));
		ans=0.5641895835*SIGN(exp(-xp*xp),x)*sum;
	}
	return ans;
}
#undef NMAX
#undef H
#undef A1
#undef A2
#undef A3
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
