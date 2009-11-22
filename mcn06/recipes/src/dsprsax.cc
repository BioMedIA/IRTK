void dsprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n)
{
	void nrerror(char error_text[]);
	unsigned long i,k;

	if (ija[1] != n+2) nrerror("dsprsax: mismatched vector and matrix");
	for (i=1;i<=n;i++) {
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<=ija[i+1]-1;k++) b[i] += sa[k]*x[ija[k]];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
