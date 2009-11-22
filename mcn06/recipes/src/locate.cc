void locate(float xx[], unsigned long n, float x, unsigned long *j)
{
	unsigned long ju,jm,jl;
	int ascnd;

	jl=0;
	ju=n+1;
	ascnd=(xx[n] > xx[1]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x > xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	*j=jl;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
