#include <math.h>

void kendl2(float **tab, int i, int j, float *tau, float *z, float *prob)
{
	float erfcc(float x);
	long nn,mm,m2,m1,lj,li,l,kj,ki,k;
	float svar,s=0.0,points,pairs,en2=0.0,en1=0.0;

	nn=i*j;
	points=tab[i][j];
	for (k=0;k<=nn-2;k++) {
		ki=(k/j);
		kj=k-j*ki;
		points += tab[ki+1][kj+1];
		for (l=k+1;l<=nn-1;l++) {
			li=l/j;
			lj=l-j*li;
			mm=(m1=li-ki)*(m2=lj-kj);
			pairs=tab[ki+1][kj+1]*tab[li+1][lj+1];
			if (mm) {
				en1 += pairs;
				en2 += pairs;
				s += (mm > 0 ? pairs : -pairs);
			} else {
				if (m1) en1 += pairs;
				if (m2) en2 += pairs;
			}
		}
	}
	*tau=s/sqrt(en1*en2);
	svar=(4.0*points+10.0)/(9.0*points*(points-1.0));
	*z=(*tau)/sqrt(svar);
	*prob=erfcc(fabs(*z)/1.4142136);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
