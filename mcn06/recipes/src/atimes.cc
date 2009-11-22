extern unsigned long ija[];
extern double sa[];

void atimes(unsigned long n, double x[], double r[], int itrnsp)
{
	void dsprsax(double sa[], unsigned long ija[], double x[], double b[],
		unsigned long n);
	void dsprstx(double sa[], unsigned long ija[], double x[], double b[],
		unsigned long n);

	if (itrnsp) dsprstx(sa,ija,x,r,n);
	else dsprsax(sa,ija,x,r,n);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
