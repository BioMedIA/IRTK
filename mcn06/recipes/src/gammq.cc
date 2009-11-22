float gammq(float a, float x)
{
	void gcf(float *gammcf, float a, float x, float *gln);
	void gser(float *gamser, float a, float x, float *gln);
	void nrerror(char error_text[]);
	float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
