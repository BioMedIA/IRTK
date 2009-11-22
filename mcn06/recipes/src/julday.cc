#include <math.h>
#define IGREG (15+31L*(10+12L*1582))

long julday(int mm, int id, int iyyy)
{
	void nrerror(char error_text[]);
	long jul;
	int ja,jy=iyyy,jm;

	if (jy == 0) nrerror("julday: there is no year zero.");
	if (jy < 0) ++jy;
	if (mm > 2) {
		jm=mm+1;
	} else {
		--jy;
		jm=mm+13;
	}
	jul = (long) (floor(365.25*jy)+floor(30.6001*jm)+id+1720995);
	if (id+31L*(mm+12L*iyyy) >= IGREG) {
		ja=(int)(0.01*jy);
		jul += 2-ja+(int) (0.25*ja);
	}
	return jul;
}
#undef IGREG
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
