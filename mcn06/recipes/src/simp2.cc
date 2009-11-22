#define EPS 1.0e-6

void simp2(float **a, int n, int l2[], int nl2, int *ip, int kp, float *q1)
{
	int k,ii,i;
	float qp,q0,q;

	*ip=0;
	for (i=1;i<=nl2;i++) {
		if (a[l2[i]+1][kp+1] < -EPS) {
			*q1 = -a[l2[i]+1][1]/a[l2[i]+1][kp+1];
			*ip=l2[i];
			for (i=i+1;i<=nl2;i++) {
				ii=l2[i];
				if (a[ii+1][kp+1] < -EPS) {
					q = -a[ii+1][1]/a[ii+1][kp+1];
					if (q < *q1) {
						*ip=ii;
						*q1=q;
					} else if (q == *q1) {
						for (k=1;k<=n;k++) {
							qp = -a[*ip+1][k+1]/a[*ip+1][kp+1];
							q0 = -a[ii+1][k+1]/a[ii+1][kp+1];
							if (q0 != qp) break;
						}
						if (q0 < qp) *ip=ii;
					}
				}
			}
		}
	}
}
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
