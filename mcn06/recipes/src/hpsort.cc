void hpsort(unsigned long n, float ra[])
{
	unsigned long i,ir,j,l;
	float rra;

	if (n < 2) return;
	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1) {
			rra=ra[--l];
		} else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				break;
			}
		}
		i=l;
		j=l+l;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) j++;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				i=j;
				j <<= 1;
			} else j=ir+1;
		}
		ra[i]=rra;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
