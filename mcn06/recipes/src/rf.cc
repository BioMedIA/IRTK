#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ERRTOL 0.08
#define TINY 1.5e-38
#define BIG 3.0e37
#define THIRD (1.0/3.0)
#define C1 (1.0/24.0)
#define C2 0.1
#define C3 (3.0/44.0)
#define C4 (1.0/14.0)

float rf(float x, float y, float z)
{
	float alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;

	if (FMIN(FMIN(x,y),z) < 0.0 || FMIN(FMIN(x+y,x+z),y+z) < TINY ||
		FMAX(FMAX(x,y),z) > BIG)
			nrerror("invalid arguments in rf");
	xt=x;
	yt=y;
	zt=z;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=THIRD*(xt+yt+zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (FMAX(FMAX(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	e2=delx*dely-delz*delz;
	e3=delx*dely*delz;
	return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
}
#undef ERRTOL
#undef TINY
#undef BIG
#undef THIRD
#undef C1
#undef C2
#undef C3
#undef C4
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
