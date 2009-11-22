void sprstp(float sa[], unsigned long ija[], float sb[], unsigned long ijb[])
{
	void iindexx(unsigned long n, long arr[], unsigned long indx[]);
	unsigned long j,jl,jm,jp,ju,k,m,n2,noff,inc,iv;
	float v;

	n2=ija[1];
	for (j=1;j<=n2-2;j++) sb[j]=sa[j];
	iindexx(ija[n2-1]-ija[1],(long *)&ija[n2-1],&ijb[n2-1]);
	jp=0;
	for (k=ija[1];k<=ija[n2-1]-1;k++) {
		m=ijb[k]+n2-1;
		sb[k]=sa[m];
		for (j=jp+1;j<=ija[m];j++) ijb[j]=k;
		jp=ija[m];
		jl=1;
		ju=n2-1;
		while (ju-jl > 1) {
			jm=(ju+jl)/2;
			if (ija[jm] > m) ju=jm; else jl=jm;
		}
		ijb[k]=jl;
	}
	for (j=jp+1;j<n2;j++) ijb[j]=ija[n2-1];
	for (j=1;j<=n2-2;j++) {
		jl=ijb[j+1]-ijb[j];
		noff=ijb[j]-1;
		inc=1;
		do {
			inc *= 3;
			inc++;
		} while (inc <= jl);
		do {
			inc /= 3;
			for (k=noff+inc+1;k<=noff+jl;k++) {
				iv=ijb[k];
				v=sb[k];
				m=k;
				while (ijb[m-inc] > iv) {
					ijb[m]=ijb[m-inc];
					sb[m]=sb[m-inc];
					m -= inc;
					if (m-noff <= inc) break;
				}
				ijb[m]=iv;
				sb[m]=v;
			}
		} while (inc > 1);
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
