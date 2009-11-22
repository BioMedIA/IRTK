void chint(float a, float b, float c[], float cint[], int n)
{
	int j;
	float sum=0.0,fac=1.0,con;

	con=0.25*(b-a);
	for (j=1;j<=n-2;j++) {
		cint[j]=con*(c[j-1]-c[j+1])/j;
		sum += fac*cint[j];
		fac = -fac;
	}
	cint[n-1]=con*c[n-2]/(n-1);
	sum += fac*cint[n-1];
	cint[0]=2.0*sum;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
