#include <math.h>
#define NRANSI
#include "nrutil.h"
#define MAXITS 200
#define EPS 1.0e-7
#define TOLF 1.0e-4
#define TOLX EPS
#define STPMX 100.0
#define TOLMIN 1.0e-6
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
	free_vector(w,1,n);free_vector(t,1,n);free_vector(s,1,n);\
	free_matrix(r,1,n,1,n);free_matrix(qt,1,n,1,n);free_vector(p,1,n);\
	free_vector(g,1,n);free_vector(fvcold,1,n);free_vector(d,1,n);\
	free_vector(c,1,n);return;}

int nn;
float *fvec;
void (*nrfuncv)(int n, float v[], float f[]);

void broydn(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []))
{
	void fdjac(int n, float x[], float fvec[], float **df,
		void (*vecfunc)(int, float [], float []));
	float fmin(float x[]);
	void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
		 float *f, float stpmax, int *check, float (*func)(float []));
	void qrdcmp(float **a, int n, float *c, float *d, int *sing);
	void qrupdt(float **r, float **qt, int n, float u[], float v[]);
	void rsolv(float **a, int n, float d[], float b[]);
	int i,its,j,k,restrt,sing,skip;
	float den,f,fold,stpmax,sum,temp,test,*c,*d,*fvcold;
	float *g,*p,**qt,**r,*s,*t,*w,*xold;

	c=vector(1,n);
	d=vector(1,n);
	fvcold=vector(1,n);
	g=vector(1,n);
	p=vector(1,n);
	qt=matrix(1,n,1,n);
	r=matrix(1,n,1,n);
	s=vector(1,n);
	t=vector(1,n);
	w=vector(1,n);
	xold=vector(1,n);
	fvec=vector(1,n);
	nn=n;
	nrfuncv=vecfunc;
	f=fmin(x);
	test=0.0;
	for (i=1;i<=n;i++)
		if (fabs(fvec[i]) > test)test=fabs(fvec[i]);
	if (test<0.01*TOLF) FREERETURN
	for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
	stpmax=STPMX*FMAX(sqrt(sum),(float)n);
	restrt=1;
	for (its=1;its<=MAXITS;its++) {
		if (restrt) {
			fdjac(n,x,fvec,r,vecfunc);
			qrdcmp(r,n,c,d,&sing);
			if (sing) nrerror("singular Jacobian in broydn");
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) qt[i][j]=0.0;
				qt[i][i]=1.0;
			}
			for (k=1;k<n;k++) {
				if (c[k]) {
					for (j=1;j<=n;j++) {
						sum=0.0;
						for (i=k;i<=n;i++)
							sum += r[i][k]*qt[i][j];
						sum /= c[k];
						for (i=k;i<=n;i++)
							qt[i][j] -= sum*r[i][k];
					}
				}
			}
			for (i=1;i<=n;i++) {
				r[i][i]=d[i];
				for (j=1;j<i;j++) r[i][j]=0.0;
			}
		} else {
			for (i=1;i<=n;i++) s[i]=x[i]-xold[i];
			for (i=1;i<=n;i++) {
				for (sum=0.0,j=i;j<=n;j++) sum += r[i][j]*s[j];
				t[i]=sum;
			}
			skip=1;
			for (i=1;i<=n;i++) {
				for (sum=0.0,j=1;j<=n;j++) sum += qt[j][i]*t[j];
				w[i]=fvec[i]-fvcold[i]-sum;
				if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i]))) skip=0;
				else w[i]=0.0;
			}
			if (!skip) {
				for (i=1;i<=n;i++) {
					for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*w[j];
					t[i]=sum;
				}
				for (den=0.0,i=1;i<=n;i++) den += SQR(s[i]);
				for (i=1;i<=n;i++) s[i] /= den;
				qrupdt(r,qt,n,t,s);
				for (i=1;i<=n;i++) {
					if (r[i][i] == 0.0) nrerror("r singular in broydn");
					d[i]=r[i][i];
				}
			}
		}
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
			g[i]=sum;
		}
		for (i=n;i>=1;i--) {
			for (sum=0.0,j=1;j<=i;j++) sum += r[j][i]*g[j];
			g[i]=sum;
		}
		for (i=1;i<=n;i++) {
			xold[i]=x[i];
			fvcold[i]=fvec[i];
		}
		fold=f;
		for (i=1;i<=n;i++) {
			for (sum=0.0,j=1;j<=n;j++) sum += qt[i][j]*fvec[j];
			p[i] = -sum;
		}
		rsolv(r,n,d,p);
		lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
		test=0.0;
		for (i=1;i<=n;i++)
			if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
		if (test < TOLF) {
			*check=0;
			FREERETURN
		}
		if (*check) {
			if (restrt) FREERETURN
			else {
				test=0.0;
				den=FMAX(f,0.5*n);
				for (i=1;i<=n;i++) {
					temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
					if (temp > test) test=temp;
				}
				if (test < TOLMIN) FREERETURN
				else restrt=1;
			}
		} else {
			restrt=0;
			test=0.0;
			for (i=1;i<=n;i++) {
				temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
				if (temp > test) test=temp;
			}
			if (test < TOLX) FREERETURN
		}
	}
	nrerror("MAXITS exceeded in broydn");
	FREERETURN
}
#undef MAXITS
#undef EPS
#undef TOLF
#undef TOLMIN
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
