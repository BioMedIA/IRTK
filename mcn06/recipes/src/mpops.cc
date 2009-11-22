#define LOBYTE(x) ((unsigned char) ((x) & 0xff))
#define HIBYTE(x) ((unsigned char) ((x) >> 8 & 0xff))


void mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int n)
{
	int j;
	unsigned short ireg=0;

	for (j=n;j>=1;j--) {
		ireg=u[j]+v[j]+HIBYTE(ireg);
		w[j+1]=LOBYTE(ireg);
	}
	w[1]=HIBYTE(ireg);
}

void mpsub(int *is, unsigned char w[], unsigned char u[], unsigned char v[],
	int n)
{
	int j;
	unsigned short ireg=256;

	for (j=n;j>=1;j--) {
		ireg=255+u[j]-v[j]+HIBYTE(ireg);
		w[j]=LOBYTE(ireg);
	}
	*is=HIBYTE(ireg)-1;
}

void mpsad(unsigned char w[], unsigned char u[], int n, int iv)
{
	int j;
	unsigned short ireg;

	ireg=256*iv;
	for (j=n;j>=1;j--) {
		ireg=u[j]+HIBYTE(ireg);
		w[j+1]=LOBYTE(ireg);
	}
	w[1]=HIBYTE(ireg);
}

void mpsmu(unsigned char w[], unsigned char u[], int n, int iv)
{
	int j;
	unsigned short ireg=0;

	for (j=n;j>=1;j--) {
		ireg=u[j]*iv+HIBYTE(ireg);
		w[j+1]=LOBYTE(ireg);
	}
	w[1]=HIBYTE(ireg);
}

void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir)
{
	int i,j;

	*ir=0;
	for (j=1;j<=n;j++) {
		i=256*(*ir)+u[j];
		w[j]=(unsigned char) (i/iv);
		*ir=i % iv;
	}
}

void mpneg(unsigned char u[], int n)
{
	int j;
	unsigned short ireg=256;

	for (j=n;j>=1;j--) {
		ireg=255-u[j]+HIBYTE(ireg);
		u[j]=LOBYTE(ireg);
	}
}

void mpmov(unsigned char u[], unsigned char v[], int n)
{
	int j;

	for (j=1;j<=n;j++) u[j]=v[j];
}

void mplsh(unsigned char u[], int n)
{
	int j;

	for (j=1;j<=n;j++) u[j]=u[j+1];
}
#undef LOBYTE
#undef HIBYTE
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
