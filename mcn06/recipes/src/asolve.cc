extern unsigned long ija[];
extern double sa[];

void asolve(unsigned long n, double b[], double x[], int itrnsp)
{
	unsigned long i;

	for(i=1;i<=n;i++) x[i]=(sa[i] != 0.0 ? b[i]/sa[i] : b[i]);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
