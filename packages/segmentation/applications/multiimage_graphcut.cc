/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkSegmentationFunction.h>

int Zcheck = 0;	
char **inputnames = NULL,*outputname = NULL;

void usage()
{
	cerr << "Usage: graphcut [NumberOfImages] [Input_1...Input_n]" << endl;
	cerr << "-atlas [n atlas1...atlasn]         atlas information for segmentation, the first one is the foreground atlas, rest is back ground atlas" << endl;
	cerr << "-numberofcomponents [n1...nn]      number of components per atlas can't be used when do not have atlas."<<endl;
	cerr << "-temporalshift [n1...nn]           temporalshift for each image"<<endl;
	cerr << "-graphcutmode [value]              Graphcut mode 1D 2D 3D 4D (1,2,3,4 with 4 connective neighbors)"<<endl;
	cerr << "-connect 0/1	                      0 neighbor connection 1 cubic connection"<<endl;
	cerr << "-dataweight [value]              Region term weight for graphcut larger region term is more important"<<endl;
	cerr << "-timeweight [value]                Time weight for boundary term in 4D graphcut case larger time is more important" <<endl;
	cerr << "-outputname [name]                 Segmentation output prefix name" << endl;
	exit(1);
}

int main( int argc, char** argv )
{
	char buffer[255];
	irtkRealImage **atlas = NULL;
	irtkImageAttributes atr;
	irtkSegmentationFunction cf;
	int i,ok,n,*c = NULL,*ts = NULL;
	int numberOfImages = 0, numberofatlas = 0, cutmode = 3,
		connect = 0;
	double x,y,z,t,dataweight = 0.5,timeweight = 1.0;

	if( argc < 3 ) usage();

	// Number Of Images
	numberOfImages = atoi(argv[1]);
	argc--;
	argv++;
	irtkRealImage **input = new irtkRealImage *[numberOfImages];
	irtkRealImage **output = new irtkRealImage *[numberOfImages];
	inputnames = new char *[numberOfImages];

	// Read image
	for (i = 0; i < numberOfImages; i++) {
		cout << argv[1] << endl;
		input[i] = new irtkRealImage(argv[1]);
		inputnames[i] = argv[1];
		argv++;
		argc--;
	}

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-atlas") == 0)) {
			argc--;
			argv++;
			numberofatlas = atoi(argv[1]);
			atlas = new irtkRealImage*[numberofatlas];
			argc--;
			argv++;
			// Read atlas for each tissue
			for (i = 0; i < numberofatlas; i++) {
				atlas[i] = new irtkRealImage;
				atlas[i]->Read(argv[1]);
				cerr << "Image " << i <<" = " << argv[1] <<endl;
				argc--;
				argv++;
			}
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-temporalshift") == 0)) {
			argc--;
			argv++;
			ts = new int[numberOfImages];
			// Read atlas for each tissue
			for (i = 0; i < numberOfImages; i++) {	
				ts[i] = atoi(argv[1]);
				cerr << "Image " << i <<" is shifted by " << argv[1] <<" frames" <<endl;
				input[i]->GetOrigin(x,y,z,t);
				input[i]->PutOrigin(x,y,z,t+ts[i]);
				argc--;
				argv++;
			}
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-numberofcomponents") == 0)) {
			argc--;
			argv++;
			c = new int[numberofatlas+1];
			// Read atlas for each tissue
			for (i = 0; i < numberofatlas+1; i++) {	
				c[i] = atoi(argv[1]);
				cerr << "Component " << i <<" has " << argv[1] <<" GMM models" <<endl;
				argc--;
				argv++;
			}
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-outputname") == 0)) {
			argc--;
			argv++;
			outputname = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-connect") == 0)) {
			argc--;
			argv++;
			connect = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-timeweight") == 0)) {
			argc--;
			argv++;
			timeweight = 1.0/atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-graphcutmode") == 0)) {
			argc--;
			argv++;
			cutmode = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-dataweight") == 0)) {
			argc--;
			argv++;
			dataweight = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	for (i = 0; i < numberOfImages; i++) {

		//find out number of time frame
		atr = input[i]->GetImageAttributes();

		// generate output file name
		if(!outputname)
			sprintf(buffer,"segtest%.2d.nii.gz",i);
		else
			sprintf(buffer,"%s%.2d.nii.gz",outputname,i);

		output[i] = new irtkRealImage(atr);

		//blur before segmentation
		cout << "Blurring for segmentation ... "; cout.flush();
		irtkGaussianBlurringWithPadding<irtkRealPixel> blurring(1, 0);
		blurring.SetInput (input[i]);
		blurring.SetOutput(input[i]);
		blurring.Run();
		cout << "done" << endl;
	}
	//evaluate threshold
	if(numberofatlas){
		cf.EvaluateGraphCut(output, input,atlas,numberOfImages,numberofatlas, timeweight,cutmode,connect,dataweight,0,c);
	}
	else{
		cerr << "multipleimage graphcut needs atlas information" << endl;
		exit(1);
	}
    for (n = 0; n < numberOfImages; n++) {
        if(!outputname)
            sprintf(buffer,"segtest%.2d.nii.gz",n);
        else
            sprintf(buffer,"%s%.2d.nii.gz",outputname,n);
        output[n]->Write(buffer);
    }

	for (i = 0; i < numberOfImages; i++) {
		delete input[i];
		delete output[i];
	}
	delete input;
	delete output;
	if(numberofatlas > 0){
		delete []atlas;
	}
}
