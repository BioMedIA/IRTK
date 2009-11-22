void chstwo(float bins1[], float bins2[], int nbins, int knstrn, float *df,
	float *chsq, float *prob)
{
	float gammq(float a, float x);
	int j;
	float temp;

	*df=nbins-knstrn;
	*chsq=0.0;
	for (j=1;j<=nbins;j++)
		if (bins1[j] == 0.0 && bins2[j] == 0.0)
			--(*df);
		else {
			temp=bins1[j]-bins2[j];
			*chsq += temp*temp/(bins1[j]+bins2[j]);
		}
	*prob=gammq(0.5*(*df),0.5*(*chsq));
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
