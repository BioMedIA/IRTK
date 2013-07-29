/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

char *input_name = NULL, *output_name = NULL;

void findminmax(double *var,int num,double &min, double &max){
	int i;
	min = 1000000;
	max = 0;
	for(i = 0; i < num; i++){
		if(var[i]<min) min = var[i];
		if(var[i]>max) max = var[i];
	}
}

void usage()
{
	cerr << "Usage: region [in] [out] <-Rx1 x1> <-Ry1 y1> <-Rz1 z1> <-Rt1 t1> "
		<< "<-Rx2 x2> <-Ry2 y2> <-Rz2 z2> <-Rt2 t2> <-landmarks points> <-ref image>" << endl
		<< "<-scale value>	value times of the specified interest region" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	bool ok;
	irtkGreyImage in, out, ref;
	irtkPointSet landmarks;
	int x1, x2, y1, t1, y2, z1, z2, t2;
	double tx[8],ty[8],tz[8],scale;
	double tx1,ty1,tz1,tx2,ty2,tz2;
	int i,refon;

	refon = 0;
	scale = 1;

	if (argc < 3) {
		usage();
	}

	input_name  = argv[1];
	argc--;
	argv++;
	output_name = argv[1];
	argc--;
	argv++;

	// Read input
	in.Read(input_name);

	// Default roi
	x1 = 0;
	y1 = 0;
	z1 = 0;
	t1 = 0;
	x2 = in.GetX();
	y2 = in.GetY();
	z2 = in.GetZ();
	t2 = in.GetT();

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
			argc--;
			argv++;
			x1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
			argc--;
			argv++;
			x2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
			argc--;
			argv++;
			y1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
			argc--;
			argv++;
			y2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
			argc--;
			argv++;
			z1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
			argc--;
			argv++;
			z2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rt1") == 0)) {
			argc--;
			argv++;
			t1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rt2") == 0)) {
			argc--;
			argv++;
			t2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-landmarks") == 0)) {
			argc--;
			argv++;
			landmarks.ReadVTK(argv[1]);
			for( i=0; i<landmarks.Size(); i++){
				in.WorldToImage(landmarks(i));
			}
			if(landmarks.Size() < 3){
				x1 = landmarks(0)._x; x2 = landmarks(1)._x;
				y1 = landmarks(0)._y; y2 = landmarks(1)._y;
				z1 = landmarks(0)._z; z2 = landmarks(1)._z;
				if (x1 > x2){
					swap(x1, x2);
				}
				if (y1 > y2){
					swap(y1, y2);
				}
				if (z1 > z2){
					swap(z1, z2);
				}
			}else{
				x1 = x2;
				x2 = 0;
				y1 = y2;
				y2 = 0;
				for( i=0; i<landmarks.Size(); i++){
					if(landmarks(i)._x < x1){
						x1 = landmarks(i)._x;
					}
					if(landmarks(i)._x > x2){
						x2 = landmarks(i)._x;
					}
					if(landmarks(i)._y < y1){
						y1 = landmarks(i)._y;
					}
					if(landmarks(i)._y > y2){
						y2 = landmarks(i)._y;
					}
				}
			}
			cout << "Using landmarks to define ROI in image coordinates. The two corners are:" << endl;
			cout << "    (" << x1 << ", " << y1 << ", " << z1 << ")" << endl;
			cout << "    (" << x2 << ", " << y2 << ", " << z2 << ")" << endl;
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && strcmp(argv[1], "-ref") == 0) {
			argc--;
			argv++;
			refon = 1;
			ok = true;
			ref.Read(argv[1]);
			//refatr = ref.GetImageAttributes();
			argc--;
			argv++;
		}
		if ((ok == false) && strcmp(argv[1], "-scale") == 0) {
			argc--;
			argv++;
			ok = true;
			scale = atof(argv[1]);
			//refatr = ref.GetImageAttributes();
			argc--;
			argv++;
		}
		if (!ok) {
			cerr << "Invalid option : " << argv[1] << endl;
			exit(1);
		}
	}

	if(refon){
		tx[0] = 0; ty[0] = 0; tz[0] = 0;
		tx[1] = 0; ty[1] = 0; tz[1] = ref.GetZ();
		tx[2] = 0; ty[2] = ref.GetY(); tz[2] = ref.GetZ();
		tx[3] = 0; ty[3] = ref.GetY(); tz[3] = 0;
		tx[4] = ref.GetX(); ty[4] = 0; tz[4] = 0;
		tx[5] = ref.GetX(); ty[5] = 0; tz[5] = ref.GetZ();
		tx[6] = ref.GetX(); ty[6] = ref.GetY(); tz[6] = 0;
		tx[7] = ref.GetX(); ty[7] = ref.GetY(); tz[7] = ref.GetZ();

		for (i=0;i<8;i++){
			ref.ImageToWorld(tx[i],ty[i],tz[i]);
			in.WorldToImage(tx[i],ty[i],tz[i]);
		}
		findminmax(tx,8,tx1,tx2);
		findminmax(ty,8,ty1,ty2);
		findminmax(tz,8,tz1,tz2);

		x1 = round(tx1); x2 = round(tx2); y1 = round(ty1);
		y2 = round(ty2); z1 = round(tz1); z2 = round(tz2);
		if(x1<0) x1 = 0; if(y1<0) y1 = 0; if(z1<0) z1 = 0;
		if(x2 == 0) x2 = 1;
		if(y2 == 0) y2 = 1;
		if(z2 == 0) z2 = 1;
		if(x2 > in.GetX()) x2 = in.GetX();
		if(y2 > in.GetY()) y2 = in.GetY();
		if(z2 > in.GetZ()) z2 = in.GetZ();
	}

	if(scale != 1.0){
		tx1 = (x1 + x2)/2.0;
		ty1 = (y1 + y2)/2.0;
		tx2 = abs(x1 - x2)/2.0;
		ty2 = abs(y1 - y2)/2.0;
		tx2 = tx2 * scale;
		ty2 = ty2 * scale;
		x1 = round(tx1 - tx2);
		x2 = round(tx1 + tx2);
		y1 = round(ty1 - ty2);
		y2 = round(ty1 + ty2);
		z1 = 0;
		z2 = in.GetZ();
		if(x1<0) x1 = 0; if(y1<0) y1 = 0; if(z1<0) z1 = 0;
		if(x2 == 0) x2 = 1;
		if(y2 == 0) y2 = 1;
		if(z2 == 0) z2 = 1;
		if(x2 > in.GetX()) x2 = in.GetX();
		if(y2 > in.GetY()) y2 = in.GetY();
		if(z2 > in.GetZ()) z2 = in.GetZ();
	}

	// Get region
	out = in.GetRegion(x1, y1, z1, t1, x2, y2, z2, t2);

	// Write region
	out.Write(output_name);

	return 0;
}
