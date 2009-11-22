#define NRANSI
#include "nrutil.h"

void wtn(float a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(float [], unsigned long, int))
{
	unsigned long i1,i2,i3,k,n,nnew,nprev=1,nt,ntot=1;
	int idim;
	float *wksp;

	for (idim=1;idim<=ndim;idim++) ntot *= nn[idim];
	wksp=vector(1,ntot);
	for (idim=1;idim<=ndim;idim++) {
		n=nn[idim];
		nnew=n*nprev;
		if (n > 4) {
			for (i2=0;i2<ntot;i2+=nnew) {
				for (i1=1;i1<=nprev;i1++) {
					for (i3=i1+i2,k=1;k<=n;k++,i3+=nprev) wksp[k]=a[i3];
					if (isign >= 0) {
						for(nt=n;nt>=4;nt >>= 1)
							(*wtstep)(wksp,nt,isign);
					} else {
						for(nt=4;nt<=n;nt <<= 1)
							(*wtstep)(wksp,nt,isign);
					}

					for (i3=i1+i2,k=1;k<=n;k++,i3+=nprev) a[i3]=wksp[k];
				}
			}
		}
		nprev=nnew;
	}
	free_vector(wksp,1,ntot);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
