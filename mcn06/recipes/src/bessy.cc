float bessy(int n, float x)
{
	float bessy0(float x);
	float bessy1(float x);
	void nrerror(char error_text[]);
	int j;
	float by,bym,byp,tox;

	if (n < 2) nrerror("Index n less than 2 in bessy");
	tox=2.0/x;
	by=bessy1(x);
	bym=bessy0(x);
	for (j=1;j<n;j++) {
		byp=j*tox*by-bym;
		bym=by;
		by=byp;
	}
	return by;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
