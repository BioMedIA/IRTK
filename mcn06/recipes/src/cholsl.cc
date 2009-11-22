void cholsl(float **a, int n, float p[], float b[], float x[])
{
	int i,k;
	float sum;

	for (i=1;i<=n;i++) {
		for (sum=b[i],k=i-1;k>=1;k--) sum -= a[i][k]*x[k];
		x[i]=sum/p[i];
	}
	for (i=n;i>=1;i--) {
		for (sum=x[i],k=i+1;k<=n;k++) sum -= a[k][i]*x[k];
		x[i]=sum/p[i];
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
