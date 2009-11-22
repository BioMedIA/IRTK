#include <stdio.h>
#define NRANSI
#include "nrutil.h"
#define IAOFF 48

void mppi(int n)
{
	void mp2dfr(unsigned char a[], unsigned char s[], int n, int *m);
	void mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int n);
	void mpinv(unsigned char u[], unsigned char v[], int n, int m);
	void mplsh(unsigned char u[], int n);
	void mpmov(unsigned char u[], unsigned char v[], int n);
	void mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
		int m);
	void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir);
	void mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int n,
		int m);
	int ir,j,m;
	unsigned char mm,*x,*y,*sx,*sxi,*t,*s,*pi;

	x=cvector(1,n+1);
	y=cvector(1,n<<1);
	sx=cvector(1,n);
	sxi=cvector(1,n);
	t=cvector(1,n<<1);
	s=cvector(1,3*n);
	pi=cvector(1,n+1);
	t[1]=2;
	for (j=2;j<=n;j++) t[j]=0;
	mpsqrt(x,x,t,n,n);
	mpadd(pi,t,x,n);
	mplsh(pi,n);
	mpsqrt(sx,sxi,x,n,n);
	mpmov(y,sx,n);
	for (;;) {
		mpadd(x,sx,sxi,n);
		mpsdv(x,&x[1],n,2,&ir);
		mpsqrt(sx,sxi,x,n,n);
		mpmul(t,y,sx,n,n);
		mpadd(&t[1],&t[1],sxi,n);
		x[1]++;
		y[1]++;
		mpinv(s,y,n,n);
		mpmul(y,&t[2],s,n,n);
		mplsh(y,n);
		mpmul(t,x,s,n,n);
		mm=t[2]-1;
		j=t[n+1]-mm;
		if (j > 1 || j < -1) {
			for (j=3;j<=n;j++) {
				if (t[j] != mm) {
					mpmul(s,pi,&t[1],n,n);
					mpmov(pi,&s[1],n);
					break;
				}
			}
			if (j <= n) continue;
		}
		printf("pi=\n");
		s[1]=pi[1]+IAOFF;
		s[2]='.';
		m=mm;
		mp2dfr(&pi[1],&s[2],n-1,&m);
		s[m+3]=0;
		printf(" %64s\n",&s[1]);
		free_cvector(pi,1,n+1);
		free_cvector(s,1,3*n);
		free_cvector(t,1,n<<1);
		free_cvector(sxi,1,n);
		free_cvector(sx,1,n);
		free_cvector(y,1,n<<1);
		free_cvector(x,1,n+1);
		return;
	}
}
#undef IAOFF
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
