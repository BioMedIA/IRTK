void hufapp(unsigned long index[], unsigned long nprob[], unsigned long n,
	unsigned long i)
{
	unsigned long j,k;

	k=index[i];
	while (i <= (n>>1)) {
		if ((j = i << 1) < n && nprob[index[j]] > nprob[index[j+1]]) j++;
		if (nprob[k] <= nprob[index[j]]) break;
		index[i]=index[j];
		i=j;
	}
	index[i]=k;
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
