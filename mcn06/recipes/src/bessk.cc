float bessk(int n, float x)
{
	float bessk0(float x);
	float bessk1(float x);
	void nrerror(char error_text[]);
	int j;
	float bk,bkm,bkp,tox;

	if (n < 2) nrerror("Index n less than 2 in bessk");
	tox=2.0/x;
	bkm=bessk0(x);
	bk=bessk1(x);
	for (j=1;j<n;j++) {
		bkp=bkm+j*tox*bk;
		bkm=bk;
		bk=bkp;
	}
	return bk;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
