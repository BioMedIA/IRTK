#include "complex.h"
#define ONE Complex(1.0,0.0)

extern fcomplex aa,bb,cc,z0,dz;

void hypdrv(float s, float yy[], float dyyds[])
{
	fcomplex z,y[3],dyds[3];

	y[1]=Complex(yy[1],yy[2]);
	y[2]=Complex(yy[3],yy[4]);
	z=Cadd(z0,RCmul(s,dz));
	dyds[1]=Cmul(y[2],dz);
	dyds[2]=Cmul(Csub(Cmul(Cmul(aa,bb),y[1]),Cmul(Csub(cc,
		Cmul(Cadd(Cadd(aa,bb),ONE),z)),y[2])),
		Cdiv(dz,Cmul(z,Csub(ONE,z))));
	dyyds[1]=dyds[1].r;
	dyyds[2]=dyds[1].i;
	dyyds[3]=dyds[2].r;
	dyyds[4]=dyds[2].i;
}
#undef ONE
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
