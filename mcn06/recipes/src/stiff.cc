#include <math.h>
#define NRANSI
#include "nrutil.h"
#define SAFETY 0.9
#define GROW 1.5
#define PGROW -0.25
#define SHRNK 0.5
#define PSHRNK (-1.0/3.0)
#define ERRCON 0.1296
#define MAXTRY 40
#define GAM (1.0/2.0)
#define A21 2.0
#define A31 (48.0/25.0)
#define A32 (6.0/25.0)
#define C21 -8.0
#define C31 (372.0/25.0)
#define C32 (12.0/5.0)
#define C41 (-112.0/125.0)
#define C42 (-54.0/125.0)
#define C43 (-2.0/5.0)
#define B1 (19.0/9.0)
#define B2 (1.0/2.0)
#define B3 (25.0/108.0)
#define B4 (125.0/108.0)
#define E1 (17.0/54.0)
#define E2 (7.0/36.0)
#define E3 0.0
#define E4 (125.0/108.0)
#define C1X (1.0/2.0)
#define C2X (-3.0/2.0)
#define C3X (121.0/50.0)
#define C4X (29.0/250.0)
#define A2X 1.0
#define A3X (3.0/5.0)

void stiff(float y[], float dydx[], int n, float *x, float htry, float eps,
	float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []))
{
	void jacobn(float x, float y[], float dfdx[], float **dfdy, int n);
	void lubksb(float **a, int n, int *indx, float b[]);
	void ludcmp(float **a, int n, int *indx, float *d);
	int i,j,jtry,*indx;
	float d,errmax,h,xsav,**a,*dfdx,**dfdy,*dysav,*err;
	float *g1,*g2,*g3,*g4,*ysav;

	indx=ivector(1,n);
	a=matrix(1,n,1,n);
	dfdx=vector(1,n);
	dfdy=matrix(1,n,1,n);
	dysav=vector(1,n);
	err=vector(1,n);
	g1=vector(1,n);
	g2=vector(1,n);
	g3=vector(1,n);
	g4=vector(1,n);
	ysav=vector(1,n);
	xsav=(*x);
	for (i=1;i<=n;i++) {
		ysav[i]=y[i];
		dysav[i]=dydx[i];
	}
	jacobn(xsav,ysav,dfdx,dfdy,n);
	h=htry;
	for (jtry=1;jtry<=MAXTRY;jtry++) {
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) a[i][j] = -dfdy[i][j];
			a[i][i] += 1.0/(GAM*h);
		}
		ludcmp(a,n,indx,&d);
		for (i=1;i<=n;i++)
			g1[i]=dysav[i]+h*C1X*dfdx[i];
		lubksb(a,n,indx,g1);
		for (i=1;i<=n;i++)
			y[i]=ysav[i]+A21*g1[i];
		*x=xsav+A2X*h;
		(*derivs)(*x,y,dydx);
		for (i=1;i<=n;i++)
			g2[i]=dydx[i]+h*C2X*dfdx[i]+C21*g1[i]/h;
		lubksb(a,n,indx,g2);
		for (i=1;i<=n;i++)
			y[i]=ysav[i]+A31*g1[i]+A32*g2[i];
		*x=xsav+A3X*h;
		(*derivs)(*x,y,dydx);
		for (i=1;i<=n;i++)
			g3[i]=dydx[i]+h*C3X*dfdx[i]+(C31*g1[i]+C32*g2[i])/h;
		lubksb(a,n,indx,g3);
		for (i=1;i<=n;i++)
			g4[i]=dydx[i]+h*C4X*dfdx[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/h;
		lubksb(a,n,indx,g4);
		for (i=1;i<=n;i++) {
			y[i]=ysav[i]+B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i];
			err[i]=E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i];
		}
		*x=xsav+h;
		if (*x == xsav) nrerror("stepsize not significant in stiff");
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(err[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) {
			*hdid=h;
#ifdef WIN32
      *hnext=(errmax > ERRCON ? SAFETY*h*pow(static_cast<double>(errmax),static_cast<double>(PGROW)) : GROW*h);
#else
			*hnext=(errmax > ERRCON ? SAFETY*h*pow(errmax,PGROW) : GROW*h);
#endif
			free_vector(ysav,1,n);
			free_vector(g4,1,n);
			free_vector(g3,1,n);
			free_vector(g2,1,n);
			free_vector(g1,1,n);
			free_vector(err,1,n);
			free_vector(dysav,1,n);
			free_matrix(dfdy,1,n,1,n);
			free_vector(dfdx,1,n);
			free_matrix(a,1,n,1,n);
			free_ivector(indx,1,n);
			return;
		} else {
#ifdef WIN32
      *hnext=SAFETY*h*pow(static_cast<double>(errmax),static_cast<double>(PSHRNK));
#else
			*hnext=SAFETY*h*pow(errmax,PSHRNK);
#endif
			h=(h >= 0.0 ? FMAX(*hnext,SHRNK*h) : FMIN(*hnext,SHRNK*h));
		}
	}
	nrerror("exceeded MAXTRY in stiff");
}
#undef SAFETY
#undef GROW
#undef PGROW
#undef SHRNK
#undef PSHRNK
#undef ERRCON
#undef MAXTRY
#undef GAM
#undef A21
#undef A31
#undef A32
#undef C21
#undef C31
#undef C32
#undef C41
#undef C42
#undef C43
#undef B1
#undef B2
#undef B3
#undef B4
#undef E1
#undef E2
#undef E3
#undef E4
#undef C1X
#undef C2X
#undef C3X
#undef C4X
#undef A2X
#undef A3X
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
