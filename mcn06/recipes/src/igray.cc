unsigned long igray(unsigned long n, int is)
{
	int ish;
	unsigned long ans,idiv;

	if (is >= 0)
		return n ^ (n >> 1);
	ish=1;
	ans=n;
	for (;;) {
		ans ^= (idiv=ans >> ish);
		if (idiv <= 1 || ish == 16) return ans;
		ish <<= 1;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
