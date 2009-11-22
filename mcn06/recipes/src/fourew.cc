#include <stdio.h>
#define SWAP(a,b) ftemp=(a);(a)=(b);(b)=ftemp

void fourew(FILE *file[5], int *na, int *nb, int *nc, int *nd)
{
	int i;
	FILE *ftemp;

	for (i=1;i<=4;i++) rewind(file[i]);
	SWAP(file[2],file[4]);
	SWAP(file[1],file[3]);
	*na=3;
	*nb=4;
	*nc=1;
	*nd=2;
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
