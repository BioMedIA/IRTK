void sprspm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
	float sc[], unsigned long ijc[])
{
	void nrerror(char error_text[]);
	unsigned long i,ijma,ijmb,j,m,ma,mb,mbb,mn;
	float sum;

	if (ija[1] != ijb[1] || ija[1] != ijc[1])
		nrerror("sprspm: sizes do not match");
	for (i=1;i<=ijc[1]-2;i++) {
		j=m=i;
		mn=ijc[i];
		sum=sa[i]*sb[i];
		for (;;) {
			mb=ijb[j];
			for (ma=ija[i];ma<=ija[i+1]-1;ma++) {
				ijma=ija[ma];
				if (ijma == j) sum += sa[ma]*sb[j];
				else {
					while (mb < ijb[j+1]) {
						ijmb=ijb[mb];
						if (ijmb == i) {
							sum += sa[i]*sb[mb++];
							continue;
						} else if (ijmb < ijma) {
							mb++;
							continue;
						} else if (ijmb == ijma) {
							sum += sa[ma]*sb[mb++];
							continue;
						}
						break;
					}
				}
			}
			for (mbb=mb;mbb<=ijb[j+1]-1;mbb++) {
				if (ijb[mbb] == i) sum += sa[i]*sb[mbb];
			}
			sc[m]=sum;
			sum=0.0;
			if (mn >= ijc[i+1]) break;
			j=ijc[m=mn++];
		}
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
