/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkImage.h>

#include <irtkResampling.h>

#ifdef HAS_CONTRIB

#include <irtkEuclideanDistanceTransform.h>
#include <irtkCityBlockDistanceTransform.h>

char *input_name = NULL, *output_name = NULL;
char *modeName = NULL;

void usage()
{
  cerr << " Usage: dmap [in] [out] <-options>"<<endl;
  cerr << " " << endl;
  cerr << " Where options are: " << endl << endl;
  cerr << "   -mode [mode]  Where [mode] can be 'euclidean' (default) or 'cityblock'." << endl << endl;

  cerr << "   -clearSlices  Set distance map to zero for completely homogeneous slices." << endl;
  cerr << "                 Used for visualisation." << endl << endl;

  cerr << " Options applying to the default Euclidean transform only are:  " << endl;
  cerr << "   <-3D/-2D>  " << endl;
  cerr << "   -radial    " << endl;
  cerr << "   -isotropic [zaxis times]" << endl << endl;
  exit(1);
}

int main(int argc, char **argv)
{
	int x, y, z, ok, radialon, isotropic, sum;
	irtkRealImage input, inputA, inputB, outputA, outputB;
	bool clearSlices = false;

	if (argc < 3) {
		usage();
	}
	radialon = 0;
	isotropic = 0;
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
		if ((ok == false) && (strcmp(argv[1], "-isotropic") == 0)) {
			argc--;
			argv++;
			if (argc > 1 && argv[1][0] != '-') {
				isotropic = atoi(argv[1]);
				argc--;
				argv++;
			} else
				isotropic = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-mode") == 0)) {
			argc--;
			argv++;
			modeName = argv[1];
			argc--;
			argv++;
			ok = true;
		}
    if ((ok == false) && (strcmp(argv[1], "-clearSlices") == 0)) {
      argc--;
      argv++;
      clearSlices = true;
      ok = true;
    }
		if (ok == false) {
			cerr << "Can't parse argument " << argv[1] << endl;
			usage();
		}
	}

	// Check if we are to do the non-default transform.
	if (modeName != NULL){
		if (strcmp(modeName, "cityblock") == 0){
			// Find City Block distance transform.

		  irtkCityBlockDistanceTransform<irtkRealPixel> cbdt;

		  // Read input
			input.Read(input_name);

			cout << "Finding city block distance transform." << endl;
			cout << "Doing inside DT" << endl;
			cbdt.SetInput(&input);
		  // Interior positive.
		  cbdt.SetOutput(&outputA);
		  cbdt.Run();

			cout << "Doing outside DT" << endl;
		  // Flip object and background voxels.
			for (z = 0; z < input.GetZ(); z++) {
				for (y = 0; y < input.GetY(); y++) {
					for (x = 0; x < input.GetX(); x++) {
						input(x, y, z) =  input(x,y,z) > 0 ? 0 : 1;
					}
				}
			}
		  cbdt.SetInput(&input);
		  // Exterior positive.
		  cbdt.SetOutput(&outputB);
		  cbdt.Run();

			for (z = 0; z < input.GetZ(); z++) {
				for (y = 0; y < input.GetY(); y++) {
					for (x = 0; x < input.GetX(); x++) {
						outputA(x, y, z) = outputB(x, y, z) - outputA(x, y, z);
					}
				}
			}
			cout << "Done." << endl;

			outputA.Write(output_name);
			exit(0);

		} else if (! (strcmp(modeName, "euclidean") == 0)  ){
			cerr << "dmap : unrecognised mode for distance transform, choices are " << endl;
			cerr << "'euclidean' or 'cityblock'. The default is euclidean if not specified." << endl;
			cerr << endl;
			usage();
			exit(1);
		}
	}

	// Default mode : Euclidean distance transform.

	if (edt == NULL) {
		edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
		(irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);
	}

	// Read input
	input.Read(input_name);

	// Threshold image
	inputA = input;
	inputB = input;
	for (z = 0; z < input.GetZ(); z++) {
		for (y = 0; y < input.GetY(); y++) {
			for (x = 0; x < input.GetX(); x++) {
				if (input(x, y, z) > 0.5) {
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


	if (clearSlices){
	  //fix the result to better visualisation
	  for (z = 0; z < input.GetZ(); z++) {
	    sum = 0;
	    for (y = 0; y < input.GetY(); y++){
	      for (x = 0; x < input.GetX(); x++){
	        sum += input.GetAsDouble(x,y,z);
	      }
	    }
	    if (sum == 0 || sum == input.GetX()*input.GetY()){
	      for (y = 0; y < input.GetY(); y++){
	        for (x = 0; x < input.GetX(); x++){
	          outputA(x, y, z) = 0;
	        }
	      }
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

	if(isotropic != 0){
		// Resample image to isotropic voxels (smalles voxel dimension)
		double xsize, ysize, zsize, size;
		outputA.GetPixelSize(&xsize, &ysize, &zsize);
		size = xsize;
		size = (size < ysize) ? size : ysize;
		size = (size < zsize) ? size : zsize;
		cerr << "Resampling image to isotropic voxel size (in mm): "
				<< size << "z axis times " << isotropic << endl;
		irtkResampling<irtkRealPixel> resampling(size, size, size*isotropic);
		irtkLinearInterpolateImageFunction interpolator;
		resampling.SetInput ((irtkRealImage*)(&outputA));
		resampling.SetOutput((irtkRealImage*)(&outputA));
		resampling.SetInterpolator(&interpolator);
		resampling.Run();
		ok = true;
		cerr << "done.."<<endl;
	}

	// Write image
	outputA.Write(output_name);

	return 0;
}

#else

int main(int argc, char **argv)
{
	cerr << "Needs to be compiled with HAS_CONTRIB" << endl;
	return 0;
}

#endif
