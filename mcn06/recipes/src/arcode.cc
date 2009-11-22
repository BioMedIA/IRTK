#include <stdio.h>
#include <stdlib.h>
#define NWK 20
#define JTRY(j,k,m) ((long)((((double)(k))*((double)(j)))/((double)(m))))

typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;

void arcode(unsigned long *ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *lcd, int isign, arithcode *acode)
{
	void arcsum(unsigned long iin[], unsigned long iout[], unsigned long ja,
		int nwk, unsigned long nrad, unsigned long nc);
	void nrerror(char error_text[]);
	int j,k;
	unsigned long ihi,ja,jh,jl,m;

	if (!isign) {
		acode->jdif=acode->nrad-1;
		for (j=NWK;j>=1;j--) {
			acode->iupb[j]=acode->nrad-1;
			acode->ilob[j]=0;
			acode->nc=j;
			if (acode->jdif > acode->minint) return;
			acode->jdif=(acode->jdif+1)*acode->nrad-1;
		}
		nrerror("NWK too small in arcode.");
	} else {
		if (isign > 0) {
			if (*ich > acode->nch) nrerror("bad ich in arcode.");
		}
		else {
			ja=(*codep)[*lcd]-acode->ilob[acode->nc];
			for (j=acode->nc+1;j<=NWK;j++) {
				ja *= acode->nrad;
				ja += ((*codep)[*lcd+j-acode->nc]-acode->ilob[j]);
			}
			ihi=acode->nch+1;
			*ich=0;
			while (ihi-(*ich) > 1) {
				m=(*ich+ihi)>>1;
				if (ja >= JTRY(acode->jdif,acode->ncumfq[m+1],acode->ncum))
					*ich=m;
				else ihi=m;
			}
			if (*ich == acode->nch) return;
		}
		jh=JTRY(acode->jdif,acode->ncumfq[*ich+2],acode->ncum);
		jl=JTRY(acode->jdif,acode->ncumfq[*ich+1],acode->ncum);
		acode->jdif=jh-jl;
		arcsum(acode->ilob,acode->iupb,jh,NWK,acode->nrad,acode->nc);
		arcsum(acode->ilob,acode->ilob,jl,NWK,acode->nrad,acode->nc);
		for (j=acode->nc;j<=NWK;j++) {
			if (*ich != acode->nch && acode->iupb[j] != acode->ilob[j]) break;
			if (*lcd > *lcode) {
				fprintf(stderr,"Reached the end of the 'code' array.\n");
				fprintf(stderr,"Attempting to expand its size.\n");
				*lcode *= 1.5;
				if ((*codep=(unsigned char *)realloc(*codep,
					(unsigned)(*lcode*sizeof(unsigned char)))) == NULL) {
					nrerror("Size expansion failed");
				}
			}
			if (isign > 0) (*codep)[*lcd]=(unsigned char)acode->ilob[j];
			++(*lcd);
		}
		if (j > NWK) return;
		acode->nc=j;
		for(j=0;acode->jdif<acode->minint;j++)
			acode->jdif *= acode->nrad;
		if (acode->nc-j < 1) nrerror("NWK too small in arcode.");
		if (j) {
			for (k=acode->nc;k<=NWK;k++) {
				acode->iupb[k-j]=acode->iupb[k];
				acode->ilob[k-j]=acode->ilob[k];
			}
		}
		acode->nc -= j;
		for (k=NWK-j+1;k<=NWK;k++) acode->iupb[k]=acode->ilob[k]=0;
	}
	return;
}
#undef NWK
#undef JTRY
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
