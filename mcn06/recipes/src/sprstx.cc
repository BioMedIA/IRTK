void sprstx(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n)
{
	void nrerror(char error_text[]);
	unsigned long i,j,k;

	if (ija[1] != n+2) nrerror("mismatched vector and matrix in sprstx");
	for (i=1;i<=n;i++) b[i]=sa[i]*x[i];
	for (i=1;i<=n;i++) {
		for (k=ija[i];k<=ija[i+1]-1;k++) {
			j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
