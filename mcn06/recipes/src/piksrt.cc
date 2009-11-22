void piksrt(int n, float arr[])
{
	int i,j;
	float a;

	for (j=2;j<=n;j++) {
		a=arr[j];
		i=j-1;
		while (i > 0 && arr[i] > a) {
			arr[i+1]=arr[i];
			i--;
		}
		arr[i+1]=a;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
