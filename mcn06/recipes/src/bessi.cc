#include <math.h>
#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

float bessi(int n, float x)
{
	float bessi0(float x);
	void nrerror(char error_text[]);
	int j;
	float bi,bim,bip,tox,ans;

	if (n < 2) nrerror("Index n less than 2 in bessi");
	if (x == 0.0)
		return 0.0;
	else {
		tox=2.0/fabs(x);
		bip=ans=0.0;
		bi=1.0;
		for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
			bim=bip+j*tox*bi;
			bip=bi;
			bi=bim;
			if (fabs(bi) > BIGNO) {
				ans *= BIGNI;
				bi *= BIGNI;
				bip *= BIGNI;
			}
			if (j == n) ans=bip;
		}
		ans *= bessi0(x)/bi;
		return x < 0.0 && (n & 1) ? -ans : ans;
	}
}
#undef ACC
#undef BIGNO
#undef BIGNI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
