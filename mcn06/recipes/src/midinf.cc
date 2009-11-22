#define FUNC(x) ((*funk)(1.0/(x))/((x)*(x)))

float midinf(float (*funk)(float), float aa, float bb, int n)
{
	float x,tnm,sum,del,ddel,b,a;
	static float s;
	int it,j;

	b=1.0/aa;
	a=1.0/bb;
	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		return (s=(s+(b-a)*sum/tnm)/3.0);
	}
}
#undef FUNC
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
