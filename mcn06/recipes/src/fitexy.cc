#include <math.h>
#define NRANSI
#include "nrutil.h"
#define POTN 1.571000
#define BIG 1.0e30
#define PI 3.14159265
#define ACC 1.0e-3

int nn;
float *xx,*yy,*sx,*sy,*ww,aa,offs;

void fitexy(float x[], float y[], int ndat, float sigx[], float sigy[],
	float *a, float *b, float *siga, float *sigb, float *chi2, float *q)
{
	void avevar(float data[], unsigned long n, float *ave, float *var);
	float brent(float ax, float bx, float cx,
		float (*f)(float), float tol, float *xmin);
	float chixy(float bang);
	void fit(float x[], float y[], int ndata, float sig[], int mwt,
		float *a, float *b, float *siga, float *sigb, float *chi2, float *q);
	float gammq(float a, float x);
	void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
		float *fc, float (*func)(float));
	float zbrent(float (*func)(float), float x1, float x2, float tol);
	int j;
	float swap,amx,amn,varx,vary,ang[7],ch[7],scale,bmn,bmx,d1,d2,r2,
		dum1,dum2,dum3,dum4,dum5;

	xx=vector(1,ndat);
	yy=vector(1,ndat);
	sx=vector(1,ndat);
	sy=vector(1,ndat);
	ww=vector(1,ndat);
	avevar(x,ndat,&dum1,&varx);
	avevar(y,ndat,&dum1,&vary);
	scale=sqrt(varx/vary);
	nn=ndat;
	for (j=1;j<=ndat;j++) {
		xx[j]=x[j];
		yy[j]=y[j]*scale;
		sx[j]=sigx[j];
		sy[j]=sigy[j]*scale;
		ww[j]=sqrt(SQR(sx[j])+SQR(sy[j]));
	}
	fit(xx,yy,nn,ww,1,&dum1,b,&dum2,&dum3,&dum4,&dum5);
	offs=ang[1]=0.0;
	ang[2]=atan(*b);
	ang[4]=0.0;
	ang[5]=ang[2];
	ang[6]=POTN;
	for (j=4;j<=6;j++) ch[j]=chixy(ang[j]);
	mnbrak(&ang[1],&ang[2],&ang[3],&ch[1],&ch[2],&ch[3],chixy);
	*chi2=brent(ang[1],ang[2],ang[3],chixy,ACC,b);
	*chi2=chixy(*b);
	*a=aa;
	*q=gammq(0.5*(nn-2),*chi2*0.5);
	for (r2=0.0,j=1;j<=nn;j++) r2 += ww[j];
	r2=1.0/r2;
	bmx=BIG;
	bmn=BIG;
	offs=(*chi2)+1.0;
	for (j=1;j<=6;j++) {
		if (ch[j] > offs) {
			d1=fabs(ang[j]-(*b));
			while (d1 >= PI) d1 -= PI;
			d2=PI-d1;
			if (ang[j] < *b) {
				swap=d1;
				d1=d2;
				d2=swap;
			}
			if (d1 < bmx) bmx=d1;
			if (d2 < bmn) bmn=d2;
		}
	}
	if (bmx < BIG) {
		bmx=zbrent(chixy,*b,*b+bmx,ACC)-(*b);
		amx=aa-(*a);
		bmn=zbrent(chixy,*b,*b-bmn,ACC)-(*b);
		amn=aa-(*a);
		*sigb=sqrt(0.5*(bmx*bmx+bmn*bmn))/(scale*SQR(cos(*b)));
		*siga=sqrt(0.5*(amx*amx+amn*amn)+r2)/scale;
	} else (*sigb)=(*siga)=BIG;
	*a /= scale;
	*b=tan(*b)/scale;
	free_vector(ww,1,ndat);
	free_vector(sy,1,ndat);
	free_vector(sx,1,ndat);
	free_vector(yy,1,ndat);
	free_vector(xx,1,ndat);
}
#undef POTN
#undef BIG
#undef PI
#undef ACC
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
