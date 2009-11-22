#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-6
#define FREEALL free_ivector(l3,1,m);free_ivector(l2,1,m);\
	free_ivector(l1,1,n+1);

void simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase,
	int izrov[], int iposv[])
{
	void simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp,
		float *bmax);
	void simp2(float **a, int n, int l2[], int nl2, int *ip, int kp, float *q1);
	void simp3(float **a, int i1, int k1, int ip, int kp);
	int i,ip,ir,is,k,kh,kp,m12,nl1,nl2;
	int *l1,*l2,*l3;
	float q1,bmax;

	if (m != (m1+m2+m3)) nrerror("Bad input constraint counts in simplx");
	l1=ivector(1,n+1);
	l2=ivector(1,m);
	l3=ivector(1,m);
	nl1=n;
	for (k=1;k<=n;k++) l1[k]=izrov[k]=k;
	nl2=m;
	for (i=1;i<=m;i++) {
		if (a[i+1][1] < 0.0) nrerror("Bad input tableau in simplx");
		l2[i]=i;
		iposv[i]=n+i;
	}
	for (i=1;i<=m2;i++) l3[i]=1;
	ir=0;
	if (m2+m3) {
		ir=1;
		for (k=1;k<=(n+1);k++) {
			q1=0.0;
			for (i=m1+1;i<=m;i++) q1 += a[i+1][k];
			a[m+2][k] = -q1;
		}
		do {
			simp1(a,m+1,l1,nl1,0,&kp,&bmax);
			if (bmax <= EPS && a[m+2][1] < -EPS) {
				*icase = -1;
				FREEALL return;
			} else if (bmax <= EPS && a[m+2][1] <= EPS) {
				m12=m1+m2+1;
				if (m12 <= m) {
					for (ip=m12;ip<=m;ip++) {
						if (iposv[ip] == (ip+n)) {
							simp1(a,ip,l1,
								nl1,1,&kp,&bmax);
							if (bmax > 0.0)
								goto one;
						}
					}
				}
				ir=0;
				--m12;
				if (m1+1 <= m12)
					for (i=m1+1;i<=m12;i++)
						if (l3[i-m1] == 1)
							for (k=1;k<=n+1;k++)
								a[i+1][k] = -a[i+1][k];
				break;
			}
			simp2(a,n,l2,nl2,&ip,kp,&q1);
			if (ip == 0) {
				*icase = -1;
				FREEALL return;
			}
	one:	simp3(a,m+1,n,ip,kp);
			if (iposv[ip] >= (n+m1+m2+1)) {
				for (k=1;k<=nl1;k++)
					if (l1[k] == kp) break;
				--nl1;
				for (is=k;is<=nl1;is++) l1[is]=l1[is+1];
				++a[m+2][kp+1];
				for (i=1;i<=m+2;i++) a[i][kp+1] = -a[i][kp+1];
			} else {
				if (iposv[ip] >= (n+m1+1)) {
					kh=iposv[ip]-m1-n;
					if (l3[kh]) {
						l3[kh]=0;
						++a[m+2][kp+1];
						for (i=1;i<=m+2;i++)
							a[i][kp+1] = -a[i][kp+1];
					}
				}
			}
			is=izrov[kp];
			izrov[kp]=iposv[ip];
			iposv[ip]=is;
		} while (ir);
	}
	for (;;) {
		simp1(a,0,l1,nl1,0,&kp,&bmax);
		if (bmax <= 0.0) {
			*icase=0;
			FREEALL return;
		}
		simp2(a,n,l2,nl2,&ip,kp,&q1);
		if (ip == 0) {
			*icase=1;
			FREEALL return;
		}
		simp3(a,m,n,ip,kp);
		is=izrov[kp];
		izrov[kp]=iposv[ip];
		iposv[ip]=is;
	}
}
#undef EPS
#undef FREEALL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
