float erffc(float x)
{
	float gammp(float a, float x);
	float gammq(float a, float x);

	return x < 0.0 ? 1.0+gammp(0.5,x*x) : gammq(0.5,x*x);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
