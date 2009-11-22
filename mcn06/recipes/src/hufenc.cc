#include <stdio.h>
#include <stdlib.h>

typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;

void hufenc(unsigned long ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *nb, huffcode *hcode)
{
	void nrerror(char error_text[]);
	int l,n;
	unsigned long k,nc;
	static unsigned long setbit[32]={0x1L,0x2L,0x4L,0x8L,0x10L,0x20L,
		0x40L,0x80L,0x100L,0x200L,0x400L,0x800L,0x1000L,0x2000L,
		0x4000L,0x8000L,0x10000L,0x20000L,0x40000L,0x80000L,0x100000L,
		0x200000L,0x400000L,0x800000L,0x1000000L,0x2000000L,0x4000000L,
		0x8000000L,0x10000000L,0x20000000L,0x40000000L,0x80000000L};

	k=ich+1;
	if (k > hcode->nch || k < 1) nrerror("ich out of range in hufenc.");
	for (n=hcode->ncod[k]-1;n>=0;n--,++(*nb)) {
		nc=(*nb >> 3);
		if (++nc >= *lcode) {
			fprintf(stderr,"Reached the end of the 'code' array.\n");
			fprintf(stderr,"Attempting to expand its size.\n");
			*lcode *= 1.5;
			if ((*codep=(unsigned char *)realloc(*codep,
				(unsigned)(*lcode*sizeof(unsigned char)))) == NULL) {
				nrerror("Size expansion failed.");
			}
		}
		l=(*nb) & 7;
		if (!l) (*codep)[nc]=0;
		if (hcode->icod[k] & setbit[n]) (*codep)[nc] |= setbit[l];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
