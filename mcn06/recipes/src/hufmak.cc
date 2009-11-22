#define NRANSI
#include "nrutil.h"

typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;

void hufmak(unsigned long nfreq[], unsigned long nchin, unsigned long *ilong,
	unsigned long *nlong, huffcode *hcode)
{
	void hufapp(unsigned long index[], unsigned long nprob[], unsigned long n,
		unsigned long i);
	int ibit;
	long node,*up;
	unsigned long j,k,*index,n,nused,*nprob;
	static unsigned long setbit[32]={0x1L,0x2L,0x4L,0x8L,0x10L,0x20L,
		0x40L,0x80L,0x100L,0x200L,0x400L,0x800L,0x1000L,0x2000L,
		0x4000L,0x8000L,0x10000L,0x20000L,0x40000L,0x80000L,0x100000L,
		0x200000L,0x400000L,0x800000L,0x1000000L,0x2000000L,0x4000000L,
		0x8000000L,0x10000000L,0x20000000L,0x40000000L,0x80000000L};

	hcode->nch=nchin;
	index=lvector(1,(long)(2*hcode->nch-1));
	up=(long *)lvector(1,(long)(2*hcode->nch-1));
	nprob=lvector(1,(long)(2*hcode->nch-1));
	for (nused=0,j=1;j<=hcode->nch;j++) {
		nprob[j]=nfreq[j];
		hcode->icod[j]=hcode->ncod[j]=0;
		if (nfreq[j]) index[++nused]=j;
	}
	for (j=nused;j>=1;j--) hufapp(index,nprob,nused,j);
	k=hcode->nch;
	while (nused > 1) {
		node=index[1];
		index[1]=index[nused--];
		hufapp(index,nprob,nused,1);
		nprob[++k]=nprob[index[1]]+nprob[node];
		hcode->left[k]=node;
		hcode->right[k]=index[1];
		up[index[1]] = -k;
		up[node]=index[1]=k;
		hufapp(index,nprob,nused,1);
	}
	up[hcode->nodemax=k]=0;
	for (j=1;j<=hcode->nch;j++) {
		if (nprob[j]) {
			for (n=0,ibit=0,node=up[j];node;node=up[node],ibit++) {
				if (node < 0) {
					n |= setbit[ibit];
					node = -node;
				}
			}
			hcode->icod[j]=n;
			hcode->ncod[j]=ibit;
		}
	}
	*nlong=0;
	for (j=1;j<=hcode->nch;j++) {
		if (hcode->ncod[j] > *nlong) {
			*nlong=hcode->ncod[j];
			*ilong=j-1;
		}
	}
	free_lvector(nprob,1,(long)(2*hcode->nch-1));
	free_lvector((unsigned long *)up,1,(long)(2*hcode->nch-1));
	free_lvector(index,1,(long)(2*hcode->nch-1));
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
