#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ERRTOL 0.05
#define TINY 1.0e-25
#define BIG 4.5e21
#define C1 (3.0/14.0)
#define C2 (1.0/6.0)
#define C3 (9.0/22.0)
#define C4 (3.0/26.0)
#define C5 (0.25*C3)
#define C6 (1.5*C4)

float rd(float x, float y, float z)
{
	float alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
		sqrtz,sum,xt,yt,zt;

	if (FMIN(x,y) < 0.0 || FMIN(x+y,z) < TINY || FMAX(FMAX(x,y),z) > BIG)
		nrerror("invalid arguments in rd");
	xt=x;
	yt=y;
	zt=z;
	sum=0.0;
	fac=1.0;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		sum += fac/(sqrtz*(zt+alamb));
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=0.2*(xt+yt+3.0*zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (FMAX(FMAX(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	ea=delx*dely;
	eb=delz*delz;
	ec=ea-eb;
	ed=ea-6.0*eb;
	ee=ed+ec+ec;
	return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
}
#undef ERRTOL
#undef TINY
#undef BIG
#undef C1
#undef C2
#undef C3
#undef C4
#undef C5
#undef C6
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
