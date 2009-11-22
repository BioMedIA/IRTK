#include <math.h>
#define ALEN(a,b,c,d) sqrt(((b)-(a))*((b)-(a))+((d)-(c))*((d)-(c)))

float revcst(float x[], float y[], int iorder[], int ncity, int n[])
{
	float xx[5],yy[5],de;
	int j,ii;

	n[3]=1 + ((n[1]+ncity-2) % ncity);
	n[4]=1 + (n[2] % ncity);
	for (j=1;j<=4;j++) {
		ii=iorder[n[j]];
		xx[j]=x[ii];
		yy[j]=y[ii];
	}
	de = -ALEN(xx[1],xx[3],yy[1],yy[3]);
	de -= ALEN(xx[2],xx[4],yy[2],yy[4]);
	de += ALEN(xx[1],xx[4],yy[1],yy[4]);
	de += ALEN(xx[2],xx[3],yy[2],yy[3]);
	return de;
}
#undef ALEN
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
