float factln(int n)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	static float a[101];

	if (n < 0) nrerror("Negative factorial in routine factln");
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else return gammln(n+1.0);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
