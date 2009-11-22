typedef struct {
	int ncof,ioff,joff;
	float *cc,*cr;
} wavefilt;

wavefilt wfilt;

void pwtset(int n)
{
	void nrerror(char error_text[]);
	int k;
	float sig = -1.0;
	static float c4[5]={0.0,0.4829629131445341,0.8365163037378079,
			0.2241438680420134,-0.1294095225512604};
	static float c12[13]={0.0,0.111540743350, 0.494623890398, 0.751133908021,
		0.315250351709,-0.226264693965,-0.129766867567,
		0.097501605587, 0.027522865530,-0.031582039318,
		0.000553842201, 0.004777257511,-0.001077301085};
	static float c20[21]={0.0,0.026670057901, 0.188176800078, 0.527201188932,
		0.688459039454, 0.281172343661,-0.249846424327,
		-0.195946274377, 0.127369340336, 0.093057364604,
		-0.071394147166,-0.029457536822, 0.033212674059,
		0.003606553567,-0.010733175483, 0.001395351747,
		0.001992405295,-0.000685856695,-0.000116466855,
		0.000093588670,-0.000013264203};
	static float c4r[5],c12r[13],c20r[21];

	wfilt.ncof=n;
	if (n == 4) {
		wfilt.cc=c4;
		wfilt.cr=c4r;
	}
	else if (n == 12) {
		wfilt.cc=c12;
		wfilt.cr=c12r;
	}
	else if (n == 20) {
		wfilt.cc=c20;
		wfilt.cr=c20r;
	}
	else nrerror("unimplemented value n in pwtset");
	for (k=1;k<=n;k++) {
		wfilt.cr[wfilt.ncof+1-k]=sig*wfilt.cc[k];
		sig = -sig;
	}
	wfilt.ioff = wfilt.joff = -(n >> 1);
}
/* (C) Copr. 1986-92 Numerical Recipes Software 5#,. */
