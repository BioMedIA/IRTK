#include "complex.h"
#define ONE Complex(1.0,0.0)

void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z, fcomplex *series,
	fcomplex *deriv)
{
	void nrerror(char error_text[]);
	int n;
	fcomplex aa,bb,cc,fac,temp;

	deriv->r=0.0;
	deriv->i=0.0;
	fac=Complex(1.0,0.0);
	temp=fac;
	aa=a;
	bb=b;
	cc=c;
	for (n=1;n<=1000;n++) {
		fac=Cmul(fac,Cmul(aa,Cdiv(bb,cc)));
		deriv->r+=fac.r;
		deriv->i+=fac.i;
		fac=Cmul(fac,RCmul(1.0/n,z));
		*series=Cadd(temp,fac);
		if (series->r == temp.r && series->i == temp.i) return;
		temp= *series;
		aa=Cadd(aa,ONE);
		bb=Cadd(bb,ONE);
		cc=Cadd(cc,ONE);

	}
	nrerror("convergence failure in hypser");
}
#undef ONE
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
