#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 20
#define TINY 1.0e-6

void qroot(float p[], int n, float *b, float *c, float eps)
{
	void poldiv(float u[], int n, float v[], int nv, float q[], float r[]);
	int iter;
	float sc,sb,s,rc,rb,r,dv,delc,delb;
	float *q,*qq,*rem;
	float d[3];

	q=vector(0,n);
	qq=vector(0,n);
	rem=vector(0,n);
	d[2]=1.0;
	for (iter=1;iter<=ITMAX;iter++) {
		d[1]=(*b);
		d[0]=(*c);
		poldiv(p,n,d,2,q,rem);
		s=rem[0];
		r=rem[1];
		poldiv(q,(n-1),d,2,qq,rem);
		sb = -(*c)*(rc = -rem[1]);
		rb = -(*b)*rc+(sc = -rem[0]);
		dv=1.0/(sb*rc-sc*rb);
		delb=(r*sc-s*rc)*dv;
		delc=(-r*sb+s*rb)*dv;
		*b += (delb=(r*sc-s*rc)*dv);
		*c += (delc=(-r*sb+s*rb)*dv);
		if ((fabs(delb) <= eps*fabs(*b) || fabs(*b) < TINY)
			&& (fabs(delc) <= eps*fabs(*c) || fabs(*c) < TINY)) {
			free_vector(rem,0,n);
			free_vector(qq,0,n);
			free_vector(q,0,n);
			return;
		}
	}
	nrerror("Too many iterations in routine qroot");
}
#undef ITMAX
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
