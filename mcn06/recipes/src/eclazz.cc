void eclazz(int nf[], int n, int (*equiv)(int, int))
{
	int kk,jj;

	nf[1]=1;
	for (jj=2;jj<=n;jj++) {
		nf[jj]=jj;
		for (kk=1;kk<=(jj-1);kk++) {
			nf[kk]=nf[nf[kk]];
			if ((*equiv)(jj,kk)) nf[nf[nf[kk]]]=jj;
		}
	}
	for (jj=1;jj<=n;jj++) nf[jj]=nf[nf[jj]];
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
