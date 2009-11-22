void fleg(float x, float pl[], int nl)
{
	int j;
	float twox,f2,f1,d;

	pl[1]=1.0;
	pl[2]=x;
	if (nl > 2) {
		twox=2.0*x;
		f2=x;
		d=1.0;
		for (j=3;j<=nl;j++) {
			f1=d++;
			f2 += twox;
			pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d;
		}
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
