#include <math.h>
#include "complex.h"
#define EPS 6.0e-8
#define EULER 0.57721566
#define MAXIT 100
#define PIBY2 1.5707963
#define FPMIN 1.0e-30
#define TMIN 2.0
#define TRUE 1
#define ONE Complex(1.0,0.0)

void cisi(float x, float *ci, float *si)
{
	void nrerror(char error_text[]);
	int i,k,odd;
	float a,err,fact,sign,sum,sumc,sums,t,term;
	fcomplex h,b,c,d,del;

	t=fabs(x);
	if (t == 0.0) {
		*si=0.0;
		*ci = -1.0/FPMIN;
		return;
	}
	if (t > TMIN) {
		b=Complex(1.0,t);
		c=Complex(1.0/FPMIN,0.0);
		d=h=Cdiv(ONE,b);
		for (i=2;i<=MAXIT;i++) {
			a = -(i-1)*(i-1);
			b=Cadd(b,Complex(2.0,0.0));
			d=Cdiv(ONE,Cadd(RCmul(a,d),b));
			c=Cadd(b,Cdiv(Complex(a,0.0),c));
			del=Cmul(c,d);
			h=Cmul(h,del);
			if (fabs(del.r-1.0)+fabs(del.i) < EPS) break;
		}
		if (i > MAXIT) nrerror("cf failed in cisi");
		h=Cmul(Complex(cos(t),-sin(t)),h);
		*ci = -h.r;
		*si=PIBY2+h.i;
	} else {
		if (t < sqrt(FPMIN)) {
			sumc=0.0;
			sums=t;
		} else {
			sum=sums=sumc=0.0;
			sign=fact=1.0;
			odd=TRUE;
			for (k=1;k<=MAXIT;k++) {
				fact *= t/k;
				term=fact/k;
				sum += sign*term;
				err=term/fabs(sum);
				if (odd) {
					sign = -sign;
					sums=sum;
					sum=sumc;
				} else {
					sumc=sum;
					sum=sums;
				}
				if (err < EPS) break;
				odd=!odd;
			}
			if (k > MAXIT) nrerror("maxits exceeded in cisi");
		}
		*si=sums;
		*ci=sumc+log(t)+EULER;
	}
	if (x < 0.0) *si = -(*si);
}
#undef EPS
#undef EULER
#undef MAXIT
#undef PIBY2
#undef FPMIN
#undef TMIN
#undef TRUE
#undef ONE
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
