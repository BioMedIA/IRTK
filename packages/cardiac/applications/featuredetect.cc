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

#include <irtkGradientImageFilter.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
	cerr << "Usage: featuredetect [in] [out]" << endl;
	cerr << "Need at least one of the following: " << endl;
	cerr << "<-intensity>            Add intensity to the output features " << endl;
	cerr << "<-gradient>             Add gradient x y z to the output features" << endl;
	cerr << "<-spatial>              Add spatial x y z to the output features" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	irtkRealImage input;
	int intensity, gradient, spatial, ok;
	intensity = false;
	gradient = false;
	spatial = false;

	if (argc < 3) {
		usage();
	}

	input_name  = argv[1];
	argc--;
	argv++;
	output_name = argv[1];
	argc--;
	argv++;

	// Parse remaining parameters
	while (argc > 1){
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-intensity") == 0)){
			argc--;
			argv++;
			intensity = true;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-gradient") == 0)){
			argc--;
			argv++;
			gradient = true;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-spatial") == 0)){
			argc--;
			argv++;
			spatial = true;
			ok = true;
		}
		if (ok == false){
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	} 

	if(intensity == false && gradient == false && spatial == false){
		usage();
	}

	// Read input
	input.Read(input_name);

	int numberoffeatures = 0;

	if(intensity == true)
		numberoffeatures++;
	if(gradient == true){
		if(input.GetZ() > 1)
			numberoffeatures += 3;
		else
			numberoffeatures += 2;
	}
	if(spatial == true){
		if(input.GetZ() > 1)
			numberoffeatures += 3;
		else
			numberoffeatures += 2;
	}

	irtkRealImage output;
	irtkImageAttributes atr = input.GetImageAttributes();
	atr._t = numberoffeatures;

	output.Initialize(atr);

	int featureindex = 0;
	if(intensity == true){
		for(int k = 0; k < input.GetZ(); k++){
			for(int j = 0; j < input.GetY(); j++){
				for(int i = 0; i < input.GetX(); i++){
					output.Put(i,j,k,featureindex,input.Get(i,j,k));
				}
			}
		}
		featureindex++;
	}
	if(gradient == true){
		irtkGradientImageFilter<short> gradient(irtkGradientImageFilter<short>::GRADIENT_VECTOR);
		irtkGreyImage tmp = input;
		irtkGreyImage *targetgradient = new irtkGreyImage();
		gradient.SetInput (&tmp);
		gradient.SetOutput(targetgradient);
		gradient.SetPadding(-1);
		gradient.Run();
		for(int k = 0; k < input.GetZ(); k++){
			for(int j = 0; j < input.GetY(); j++){
				for(int i = 0; i < input.GetX(); i++){
					output.Put(i,j,k,featureindex,targetgradient->GetAsDouble(i,j,k,0));
					output.Put(i,j,k,featureindex+1,targetgradient->GetAsDouble(i,j,k,1));
					if(input.GetZ() > 1)
						output.Put(i,j,k,featureindex+2,targetgradient->GetAsDouble(i,j,k,2));
				}
			}
		}
		delete targetgradient;
		if(input.GetZ() > 1)
			featureindex += 3;
		else
			featureindex += 2;
	}
	if(spatial == true){
		for(int k = 0; k < input.GetZ(); k++){
			for(int j = 0; j < input.GetY(); j++){
				for(int i = 0; i < input.GetX(); i++){
					double x,y,z;
					x = i;
					y = j;
					z = k;
					input.ImageToWorld(x,y,z);
					output.Put(i,j,k,featureindex,x);
					output.Put(i,j,k,featureindex+1,y);
					if(input.GetZ() > 1)
						output.Put(i,j,k,featureindex+2,z);
				}
			}
		}
		if(input.GetZ() > 1)
			featureindex += 3;
		else
			featureindex += 2;
	}

	// Write image
	output.Write(output_name);

	return 0;
}
