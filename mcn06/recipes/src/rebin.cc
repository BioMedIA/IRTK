void rebin(float rc, int nd, float r[], float xin[], float xi[])
{
	int i,k=0;
	float dr=0.0,xn=0.0,xo;

	for (i=1;i<nd;i++) {
		while (rc > dr) {
			dr += r[++k];
			xo=xn;
			xn=xi[k];
		}
		dr -= rc;
		xin[i]=xn-(xn-xo)*dr/r[k];
	}
	for (i=1;i<nd;i++) xi[i]=xin[i];
	xi[nd]=1.0;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
