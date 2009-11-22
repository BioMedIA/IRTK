#include <math.h>
#define IGREG 2299161

void caldat(long julian, int *mm, int *id, int *iyyy)
{
	long ja,jalpha,jb,jc,jd,je;

	if (julian >= IGREG) {
		jalpha=(long)(((float) (julian-1867216)-0.25)/36524.25);
		ja=julian+1+jalpha-(long) (0.25*jalpha);
	} else
		ja=julian;
	jb=ja+1524;
	jc=(long)(6680.0+((float) (jb-2439870)-122.1)/365.25);
	jd=(long)(365*jc+(0.25*jc));
	je=(long)((jb-jd)/30.6001);
	*id=jb-jd-(long) (30.6001*je);
	*mm=je-1;
	if (*mm > 12) *mm -= 12;
	*iyyy=jc-4715;
	if (*mm > 2) --(*iyyy);
	if (*iyyy <= 0) --(*iyyy);
}
#undef IGREG
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
