#include <stdio.h>
#define ISCR 60
#define JSCR 21
#define BLANK ' '
#define ZERO '-'
#define YY 'l'
#define XX '-'
#define FF 'x'

void scrsho(float (*fx)(float))
{
	int jz,j,i;
	float ysml,ybig,x2,x1,x,dyj,dx,y[ISCR+1];
	char scr[ISCR+1][JSCR+1];

	for (;;) {
		printf("\nEnter x1 x2 (x1=x2 to stop):\n");
		scanf("%f %f",&x1,&x2);
		if (x1 == x2) break;
		for (j=1;j<=JSCR;j++)
			scr[1][j]=scr[ISCR][j]=YY;
		for (i=2;i<=(ISCR-1);i++) {
			scr[i][1]=scr[i][JSCR]=XX;
			for (j=2;j<=(JSCR-1);j++)
				scr[i][j]=BLANK;
		}
		dx=(x2-x1)/(ISCR-1);
		x=x1;
		ysml=ybig=0.0;
		for (i=1;i<=ISCR;i++) {
			y[i]=(*fx)(x);
			if (y[i] < ysml) ysml=y[i];
			if (y[i] > ybig) ybig=y[i];
			x += dx;
		}
		if (ybig == ysml) ybig=ysml+1.0;
		dyj=(JSCR-1)/(ybig-ysml);
		jz=1-(int) (ysml*dyj);
		for (i=1;i<=ISCR;i++) {
			scr[i][jz]=ZERO;
			j=1+(int) ((y[i]-ysml)*dyj);
			scr[i][j]=FF;
		}
		printf(" %10.3f ",ybig);
		for (i=1;i<=ISCR;i++) printf("%c",scr[i][JSCR]);
		printf("\n");
		for (j=(JSCR-1);j>=2;j--) {
			printf("%12s"," ");
			for (i=1;i<=ISCR;i++) printf("%c",scr[i][j]);
			printf("\n");
		}
		printf(" %10.3f ",ysml);
		for (i=1;i<=ISCR;i++) printf("%c",scr[i][1]);
		printf("\n");
		printf("%8s %10.3f %44s %10.3f\n"," ",x1," ",x2);
	}
}
#undef ISCR
#undef JSCR
#undef BLANK
#undef ZERO
#undef YY
#undef XX
#undef FF
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
