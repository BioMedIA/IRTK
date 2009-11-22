void chder(float a, float b, float c[], float cder[], int n)
{
	int j;
	float con;

	cder[n-1]=0.0;
	cder[n-2]=2*(n-1)*c[n-1];
	for (j=n-3;j>=0;j--)
		cder[j]=cder[j+2]+2*(j+1)*c[j+1];
	con=2.0/(b-a);
	for (j=0;j<n;j++)
		cder[j] *= con;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
