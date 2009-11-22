#include <math.h>

void sprstm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
	float thresh, unsigned long nmax, float sc[], unsigned long ijc[])
{
	void nrerror(char error_text[]);
	unsigned long i,ijma,ijmb,j,k,ma,mb,mbb;
	float sum;

	if (ija[1] != ijb[1]) nrerror("sprstm: sizes do not match");
	ijc[1]=k=ija[1];
	for (i=1;i<=ija[1]-2;i++) {
		for (j=1;j<=ijb[1]-2;j++) {
			if (i == j) sum=sa[i]*sb[j]; else sum=0.0e0;
			mb=ijb[j];
			for (ma=ija[i];ma<=ija[i+1]-1;ma++) {
				ijma=ija[ma];
				if (ijma == j) sum += sa[ma]*sb[j];
				else {
					while (mb < ijb[j+1]) {
						ijmb=ijb[mb];
						if (ijmb == i) {
							sum += sa[i]*sb[mb++];
							continue;
						} else if (ijmb < ijma) {
							mb++;
							continue;
						} else if (ijmb == ijma) {
							sum += sa[ma]*sb[mb++];
							continue;
						}
						break;
					}
				}
			}
			for (mbb=mb;mbb<=ijb[j+1]-1;mbb++) {
				if (ijb[mbb] == i) sum += sa[i]*sb[mbb];
			}
			if (i == j) sc[i]=sum;
			else if (fabs(sum) > thresh) {
				if (k > nmax) nrerror("sprstm: nmax too small");
				sc[k]=sum;
				ijc[k++]=j;
			}
		}
		ijc[i+1]=k;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
