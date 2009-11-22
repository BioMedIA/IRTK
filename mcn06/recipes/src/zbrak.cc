void zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[],
	float xb2[], int *nb)
{
	int nbb,i;
	float x,fp,fc,dx;

	nbb=0;
	dx=(x2-x1)/n;
	fp=(*fx)(x=x1);
	for (i=1;i<=n;i++) {
		fc=(*fx)(x += dx);
		if (fc*fp < 0.0) {
			xb1[++nbb]=x-dx;
			xb2[nbb]=x;
			if(*nb == nbb) return;

		}
		fp=fc;
	}
	*nb = nbb;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
