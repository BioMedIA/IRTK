void ddpoly(float c[], int nc, float x, float pd[], int nd)
{
	int nnd,j,i;
	float cnst=1.0;

	pd[0]=c[nc];
	for (j=1;j<=nd;j++) pd[j]=0.0;
	for (i=nc-1;i>=0;i--) {
		nnd=(nd < (nc-i) ? nd : nc-i);
		for (j=nnd;j>=1;j--)
			pd[j]=pd[j]*x+pd[j-1];
		pd[0]=pd[0]*x+c[i];
	}
	for (i=2;i<=nd;i++) {
		cnst *= i;
		pd[i] *= cnst;
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
