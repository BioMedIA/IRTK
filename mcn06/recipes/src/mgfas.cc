#define NRANSI
#include "nrutil.h"
#define NPRE 1
#define NPOST 1
#define ALPHA 0.33
#define NGMAX 15

void mgfas(double **u, int n, int maxcyc)
{
	double anorm2(double **a, int n);
	void copy(double **aout, double **ain, int n);
	void interp(double **uf, double **uc, int nf);
	void lop(double **out, double **u, int n);
	void matadd(double **a, double **b, double **c, int n);
	void matsub(double **a, double **b, double **c, int n);
	void relax2(double **u, double **rhs, int n);
	void rstrct(double **uc, double **uf, int nc);
	void slvsm2(double **u, double **rhs);
	unsigned int j,jcycle,jj,jm1,jpost,jpre,nf,ng=0,ngrid,nn;
	double **irho[NGMAX+1],**irhs[NGMAX+1],**itau[NGMAX+1],
		**itemp[NGMAX+1],**iu[NGMAX+1];
	double res,trerr;

	nn=n;
	while (nn >>= 1) ng++;
	if (n != 1+(1L << ng)) nrerror("n-1 must be a power of 2 in mgfas.");
	if (ng > NGMAX) nrerror("increase NGMAX in mglin.");
	nn=n/2+1;
	ngrid=ng-1;
	irho[ngrid]=dmatrix(1,nn,1,nn);
	rstrct(irho[ngrid],u,nn);
	while (nn > 3) {
		nn=nn/2+1;
		irho[--ngrid]=dmatrix(1,nn,1,nn);
		rstrct(irho[ngrid],irho[ngrid+1],nn);
	}
	nn=3;
	iu[1]=dmatrix(1,nn,1,nn);
	irhs[1]=dmatrix(1,nn,1,nn);
	itau[1]=dmatrix(1,nn,1,nn);
	itemp[1]=dmatrix(1,nn,1,nn);
	slvsm2(iu[1],irho[1]);
	free_dmatrix(irho[1],1,nn,1,nn);
	ngrid=ng;
	for (j=2;j<=ngrid;j++) {
		nn=2*nn-1;
		iu[j]=dmatrix(1,nn,1,nn);
		irhs[j]=dmatrix(1,nn,1,nn);
		itau[j]=dmatrix(1,nn,1,nn);
		itemp[j]=dmatrix(1,nn,1,nn);
		interp(iu[j],iu[j-1],nn);
		copy(irhs[j],(j != ngrid ? irho[j] : u),nn);
		for (jcycle=1;jcycle<=maxcyc;jcycle++) {
		nf=nn;
			for (jj=j;jj>=2;jj--) {
				for (jpre=1;jpre<=NPRE;jpre++)
					relax2(iu[jj],irhs[jj],nf);
				lop(itemp[jj],iu[jj],nf);
				nf=nf/2+1;
				jm1=jj-1;
				rstrct(itemp[jm1],itemp[jj],nf);
				rstrct(iu[jm1],iu[jj],nf);
				lop(itau[jm1],iu[jm1],nf);
				matsub(itau[jm1],itemp[jm1],itau[jm1],nf);
				if (jj == j)
					trerr=ALPHA*anorm2(itau[jm1],nf);
				rstrct(irhs[jm1],irhs[jj],nf);
				matadd(irhs[jm1],itau[jm1],irhs[jm1],nf);
			}
			slvsm2(iu[1],irhs[1]);
			nf=3;
			for (jj=2;jj<=j;jj++) {
			jm1=jj-1;
			rstrct(itemp[jm1],iu[jj],nf);
			matsub(iu[jm1],itemp[jm1],itemp[jm1],nf);
			nf=2*nf-1;
			interp(itau[jj],itemp[jm1],nf);
			matadd(iu[jj],itau[jj],iu[jj],nf);
			for (jpost=1;jpost<=NPOST;jpost++)
				relax2(iu[jj],irhs[jj],nf);
			}
			lop(itemp[j],iu[j],nf);
			matsub(itemp[j],irhs[j],itemp[j],nf);
			res=anorm2(itemp[j],nf);
			if (res < trerr) break;
		}
	}
	copy(u,iu[ngrid],n);
	for (nn=n,j=ng;j>=1;j--,nn=nn/2+1) {
		free_dmatrix(itemp[j],1,nn,1,nn);
		free_dmatrix(itau[j],1,nn,1,nn);
		free_dmatrix(irhs[j],1,nn,1,nn);
		free_dmatrix(iu[j],1,nn,1,nn);
		if (j != ng && j != 1) free_dmatrix(irho[j],1,nn,1,nn);
	}
}
#undef NGMAX
#undef NPRE
#undef NPOST
#undef ALPHA
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
