float ran4(long *idum)
{
	void psdes(unsigned long *lword, unsigned long *irword);
	unsigned long irword,itemp,lword;
	static long idums = 0;
#if defined(vax) || defined(_vax_) || defined(__vax__) || defined(VAX)
	static unsigned long jflone = 0x00004080;
	static unsigned long jflmsk = 0xffff007f;
#else
	static unsigned long jflone = 0x3f800000;
	static unsigned long jflmsk = 0x007fffff;
#endif

	if (*idum < 0) {
		idums = -(*idum);
		*idum=1;
	}
	irword=(*idum);
	lword=idums;
	psdes(&lword,&irword);
	itemp=jflone | (jflmsk & irword);
	++(*idum);
	return (*(float *)&itemp)-1.0;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
