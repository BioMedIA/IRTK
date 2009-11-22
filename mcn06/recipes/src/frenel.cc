#include <math.h>
#include "complex.h"
#define EPS 6.0e-8
#define MAXIT 100
#define FPMIN 1.0e-30
#define XMIN 1.5
#define PI 3.1415927
#define PIBY2 (PI/2.0)
#define TRUE 1
#define ONE Complex(1.0,0.0)

void frenel(float x, float *s, float *c)
{
	void nrerror(char error_text[]);
	int k,n,odd;
	float a,ax,fact,pix2,sign,sum,sumc,sums,term,test;
	fcomplex b,cc,d,h,del,cs;

	ax=fabs(x);
	if (ax < sqrt(FPMIN)) {
		*s=0.0;
		*c=ax;
	} else if (ax <= XMIN) {
		sum=sums=0.0;
		sumc=ax;
		sign=1.0;
		fact=PIBY2*ax*ax;
		odd=TRUE;
		term=ax;
		n=3;
		for (k=1;k<=MAXIT;k++) {
			term *= fact/k;
			sum += sign*term/n;
			test=fabs(sum)*EPS;
			if (odd) {
				sign = -sign;
				sums=sum;
				sum=sumc;
			} else {
				sumc=sum;
				sum=sums;
			}
			if (term < test) break;
			odd=!odd;
			n += 2;
		}
		if (k > MAXIT) nrerror("series failed in frenel");
		*s=sums;
		*c=sumc;
	} else {
		pix2=PI*ax*ax;
		b=Complex(1.0,-pix2);
		cc=Complex(1.0/FPMIN,0.0);
		d=h=Cdiv(ONE,b);
		n = -1;
		for (k=2;k<=MAXIT;k++) {
			n += 2;
			a = -n*(n+1);
			b=Cadd(b,Complex(4.0,0.0));
			d=Cdiv(ONE,Cadd(RCmul(a,d),b));
			cc=Cadd(b,Cdiv(Complex(a,0.0),cc));
			del=Cmul(cc,d);
			h=Cmul(h,del);
			if (fabs(del.r-1.0)+fabs(del.i) < EPS) break;
		}
		if (k > MAXIT) nrerror("cf failed in frenel");
		h=Cmul(Complex(ax,-ax),h);
		cs=Cmul(Complex(0.5,0.5),
			Csub(ONE,Cmul(Complex(cos(0.5*pix2),sin(0.5*pix2)),h)));
		*c=cs.r;
		*s=cs.i;
	}
	if (x < 0.0) {
		*c = -(*c);
		*s = -(*s);
	}
}
#undef EPS
#undef MAXIT
#undef FPMIN
#undef XMIN
#undef PI
#undef PIBY2
#undef TRUE
#undef ONE
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
