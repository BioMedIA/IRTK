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

#include <irtkResampling.h>

#ifdef HAS_CONTRIB

#include <irtkEuclideanDistanceTransform.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
	cerr << " Usage: label2distancemap [in] [out] <-options>"<<endl;
	cerr << "   -nobackground    " << endl;
	cerr << "   <-3D/-2D>        " << endl;
	cerr << "   -radial          " << endl;
	cerr << "   -tradial         " << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int x, y, z, ok, radialon, sum, background;
	irtkRealImage input, inputA, inputB, outputA, outputB, output;

	if (argc < 3) {
		usage();
	}
	background = true;
	radialon = 0;
	sum = 0;

	input_name  = argv[1];
	argc--;
	argv++;
	output_name = argv[1];
	argc--;
	argv++;

	// Default mode
	irtkEuclideanDistanceTransform<irtkRealPixel> *edt = NULL;
	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-2D") == 0)) {
			argc--;
			argv++;
			edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
				(irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform2D);
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-nobackground") == 0)) {
			argc--;
			argv++;
			background = false;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-3D") == 0)) {
			argc--;
			argv++;
			edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
				(irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-radial") == 0)) {
			argc--;
			argv++;
			radialon = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-tradial") == 0)) {
			argc--;
			argv++;
			radialon = 2;
			ok = true;
		}
		if (ok == false) {
			cerr << "Can't parse argument " << argv[1] << endl;
			usage();
		}
	}

	// Default mode : Euclidean distance transform.

	if (edt == NULL) {
		edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
			(irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);
	}

	// Read input
	input.Read(input_name);

	double min,max;
	int count,labelcount = 0;

	input.GetMinMax(&min,&max);

	if(background == false){
		min = min + 1;
	}

	for(int i = min; i <= max; i++){
		count = 0;
		for (z = 0; z < input.GetZ(); z++) {
			for (y = 0; y < input.GetY(); y++) {
				for (x = 0; x < input.GetX(); x++) {
					if (input(x, y, z) == i) {
						count++;
					}				
				}
			}
		}
		if(count > 0){
			labelcount++;
		}
	}

	irtkImageAttributes atr = input.GetImageAttributes();
	atr._t = labelcount;
	output.Initialize(atr);

	labelcount = 0;
	for(int i = min; i <= max; i++){
		count = 0;
		for (z = 0; z < input.GetZ(); z++) {
			for (y = 0; y < input.GetY(); y++) {
				for (x = 0; x < input.GetX(); x++) {
					if (input(x, y, z) == i) {
						count++;
					}				
				}
			}
		}
		if(count > 0){

			// Threshold image
			inputA = input;
			inputB = input;
			for (z = 0; z < input.GetZ(); z++) {
				for (y = 0; y < input.GetY(); y++) {
					for (x = 0; x < input.GetX(); x++) {
						if (input(x, y, z) == i) {
							inputA(x, y, z) = 1;
							inputB(x, y, z) = 0;
						} else {
							inputA(x, y, z) = 0;
							inputB(x, y, z) = 1;
						}
					}
				}
			}

			// Calculate EDT
			cout << "Finding Euclidean distance transform." << endl;
			cout << "Doing outside DT" << endl;
			edt->SetInput (& inputA);
			edt->SetOutput(&outputA);
			edt->Run();
			cout << "Doing inside DT" << endl;
			edt->SetInput (& inputB);
			edt->SetOutput(&outputB);
			edt->Run();

			for (z = 0; z < input.GetZ(); z++) {
				for (y = 0; y < input.GetY(); y++) {
					for (x = 0; x < input.GetX(); x++) {
						outputA(x, y, z)  = sqrt(outputA(x, y, z)) - sqrt(outputB(x, y, z));
					}
				}
			}

			if(radialon == 1){
				edt->SetInput(& outputA);
				edt->SetOutput(& outputA);
				edt->Radial();
			}else if(radialon == 2){
				edt->SetInput(& outputA);
				edt->SetOutput(& outputA);
				edt->TRadial();
			}

			// Write image
			for (z = 0; z < input.GetZ(); z++) {
				for (y = 0; y < input.GetY(); y++) {
					for (x = 0; x < input.GetX(); x++) {
						output.PutAsDouble(x,y,z,labelcount,outputA.GetAsDouble(x,y,z));
					}
				}
			}

			labelcount++;
		}
	}

	output.Write(output_name);

	return 0;
}

#else

int main(int argc, char **argv)
{
	cerr << "Needs to be compiled with HAS_CONTRIB" << endl;
	return 0;
}

#endif
