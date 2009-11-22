#define NRANSI
#include "nrutil.h"
#define NPRE 1
#define NPOST 1
#define NGMAX 15

void mglin(double **u, int n, int ncycle)
{
	void addint(double **uf, double **uc, double **res, int nf);
	void copy(double **aout, double **ain, int n);
	void fill0(double **u, int n);
	void interp(double **uf, double **uc, int nf);
	void relax(double **u, double **rhs, int n);
	void resid(double **res, double **u, double **rhs, int n);
	void rstrct(double **uc, double **uf, int nc);
	void slvsml(double **u, double **rhs);
	unsigned int j,jcycle,jj,jpost,jpre,nf,ng=0,ngrid,nn;
	double **ires[NGMAX+1],**irho[NGMAX+1],**irhs[NGMAX+1],**iu[NGMAX+1];

	nn=n;
	while (nn >>= 1) ng++;
	if (n != 1+(1L << ng)) nrerror("n-1 must be a power of 2 in mglin.");
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
	slvsml(iu[1],irho[1]);
	free_dmatrix(irho[1],1,nn,1,nn);
	ngrid=ng;
	for (j=2;j<=ngrid;j++) {
		nn=2*nn-1;
		iu[j]=dmatrix(1,nn,1,nn);
		irhs[j]=dmatrix(1,nn,1,nn);
		ires[j]=dmatrix(1,nn,1,nn);
		interp(iu[j],iu[j-1],nn);
		copy(irhs[j],(j != ngrid ? irho[j] : u),nn);
		for (jcycle=1;jcycle<=ncycle;jcycle++) {
			nf=nn;
			for (jj=j;jj>=2;jj--) {
			for (jpre=1;jpre<=NPRE;jpre++)
				relax(iu[jj],irhs[jj],nf);
			resid(ires[jj],iu[jj],irhs[jj],nf);
			nf=nf/2+1;
			rstrct(irhs[jj-1],ires[jj],nf);
			fill0(iu[jj-1],nf);
			}
			slvsml(iu[1],irhs[1]);
			nf=3;
			for (jj=2;jj<=j;jj++) {
			nf=2*nf-1;
			addint(iu[jj],iu[jj-1],ires[jj],nf);
			for (jpost=1;jpost<=NPOST;jpost++)
				relax(iu[jj],irhs[jj],nf);
			}
		}
	}
	copy(u,iu[ngrid],n);
	for (nn=n,j=ng;j>=2;j--,nn=nn/2+1) {
		free_dmatrix(ires[j],1,nn,1,nn);
		free_dmatrix(irhs[j],1,nn,1,nn);
		free_dmatrix(iu[j],1,nn,1,nn);
		if (j != ng) free_dmatrix(irho[j],1,nn,1,nn);
	}
	free_dmatrix(irhs[1],1,3,1,3);
	free_dmatrix(iu[1],1,3,1,3);
}
#undef NPRE
#undef NPOST
#undef NGMAX
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
