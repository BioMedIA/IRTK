#define NRANSI
#include "nrutil.h"
#define M 64
#define BIG 1.0e30
#define FREEALL free_vector(sel,1,M+2);free_lvector(isel,1,M+2);

float selip(unsigned long k, unsigned long n, float arr[])
{
	void shell(unsigned long n, float a[]);
	unsigned long i,j,jl,jm,ju,kk,mm,nlo,nxtmm,*isel;
	float ahi,alo,sum,*sel;

	if (k < 1 || k > n || n <= 0) nrerror("bad input to selip");
	isel=lvector(1,M+2);
	sel=vector(1,M+2);
	kk=k;
	ahi=BIG;
	alo = -BIG;
	for (;;) {
		mm=nlo=0;
		sum=0.0;
		nxtmm=M+1;
		for (i=1;i<=n;i++) {
			if (arr[i] >= alo && arr[i] <= ahi) {
				mm++;
				if (arr[i] == alo) nlo++;
				if (mm <= M) sel[mm]=arr[i];
				else if (mm == nxtmm) {
					nxtmm=mm+mm/M;
					sel[1 + ((i+mm+kk) % M)]=arr[i];
				}
				sum += arr[i];
			}
		}
		if (kk <= nlo) {
			FREEALL
			return alo;
		}
		else if (mm <= M) {
			shell(mm,sel);
			ahi = sel[kk];
			FREEALL
			return ahi;
		}
		sel[M+1]=sum/mm;
		shell(M+1,sel);
		sel[M+2]=ahi;
		for (j=1;j<=M+2;j++) isel[j]=0;
		for (i=1;i<=n;i++) {
			if (arr[i] >= alo && arr[i] <= ahi) {
				jl=0;
				ju=M+2;
				while (ju-jl > 1) {
					jm=(ju+jl)/2;
					if (arr[i] >= sel[jm]) jl=jm;
					else ju=jm;
				}
				isel[ju]++;
			}
		}
		j=1;
		while (kk > isel[j]) {
			alo=sel[j];
			kk -= isel[j++];
		}
		ahi=sel[j];
	}
}
#undef M
#undef BIG
#undef FREEALL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
