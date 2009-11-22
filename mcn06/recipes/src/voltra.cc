#define NRANSI
#include "nrutil.h"

void voltra(int n, int m, float t0, float h, float *t, float **f,
	float (*g)(int, float), float (*ak)(int, int, float, float))
{
	void lubksb(float **a, int n, int *indx, float b[]);
	void ludcmp(float **a, int n, int *indx, float *d);
	int i,j,k,l,*indx;
	float d,sum,**a,*b;

	indx=ivector(1,m);
	a=matrix(1,m,1,m);
	b=vector(1,m);
	t[1]=t0;
	for (k=1;k<=m;k++) f[k][1]=(*g)(k,t[1]);
	for (i=2;i<=n;i++) {
		t[i]=t[i-1]+h;
		for (k=1;k<=m;k++) {
			sum=(*g)(k,t[i]);
			for (l=1;l<=m;l++) {
				sum += 0.5*h*(*ak)(k,l,t[i],t[1])*f[l][1];
				for (j=2;j<i;j++)
					sum += h*(*ak)(k,l,t[i],t[j])*f[l][j];
				a[k][l]=(k == l)-0.5*h*(*ak)(k,l,t[i],t[i]);
			}
			b[k]=sum;
		}
		ludcmp(a,m,indx,&d);
		lubksb(a,m,indx,b);
		for (k=1;k<=m;k++) f[k][i]=b[k];
	}
	free_vector(b,1,m);
	free_matrix(a,1,m,1,m);
	free_ivector(indx,1,m);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
