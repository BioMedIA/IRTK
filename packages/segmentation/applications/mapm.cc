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

#include <irtkSegmentationFunction.h>

void usage()
{
	cerr << "Usage: mapm [N] [atlasfile] [target] [patchradius] [outputimage]" << endl;
	cerr << "-searchradius         [0-1] random search radius in the image, 1 means whole image" << endl;
	cerr << "-nnfiterations        [N] number of iterations of the multi-atlas patchmatch" << endl;
	cerr << "-debug                open debug mode, output intermedia results" << endl;
	exit(1);

}
int main(int argc, char **argv){

	if (argc < 6) {
		usage();
	}

	int nAtlases = atoi(argv[1]);
	cout << "N Atlases: " << nAtlases << endl;
	argv++;
	argc--;

	irtkGreyImage ** atlases = NULL;
	char *output_name = NULL;
	double xsize,ysize,zsize;
	int em_iterations = 1;
	int nnf_iterations = 40;
	double searchradius = 0.1;
	int debug = false;

	atlases = new irtkGreyImage*[nAtlases];

	ifstream infile;
	string line = "";
	infile.open (argv[1]);

	argv++;
	argc--;

	for(int i = 0; i < nAtlases; i++){
		if(!infile.eof()){
			getline(infile,line);
			cout << "Reading atlas " << line << endl;
			atlases[i] = new irtkGreyImage((char*)line.c_str());
		}else{
			cout << "Not enough atlases, should be " << nAtlases << " actually " << i+1 << endl;
			break;
		}
	}

	infile.close();

	irtkGreyImage image;
	cout << "Reading image " << argv[1] << endl;
	image.Read(argv[1]);
	argc--;
	argv++;
	int patchSize = atoi(argv[1]);
	argc--;
	argv++;
	cout << "Patch radius: " << patchSize << endl;

	char * outImageName = argv[1];
	argc--;
	argv++;

	char * outLabelName = NULL;

	while (argc > 1){
		int ok = false;
		if ((ok == false) && (strcmp(argv[1], "-output") == 0)){
			argv++;
			argc--;
			output_name = argv[1];
			argv++;
			argc--;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-nnfiterations") == 0)){
			argv++;
			argc--;
			nnf_iterations = atoi(argv[1]);
			argv++;
			argc--;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-searchradius") == 0)){
			argv++;
			argc--;
			searchradius = atof(argv[1]);
			argv++;
			argc--;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-debug") == 0)){
			argv++;
			argc--;
			debug = true;
			ok = true;
		}
		if (ok == false){
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	cout << "Creating patchmatch..."<<endl;

	irtkMAPatchMatch *patchmatch = new irtkMAPatchMatch(&image, atlases, patchSize, nAtlases, 1);

	patchmatch->setDebug(debug);

	cout << "Creating patchmatch done"<<endl;

	cout << "Start optimization..."<<endl;

	patchmatch->setRandomrate(searchradius);

	patchmatch->run(nnf_iterations);

	cout << "Optimization done"<<endl;
		
	patchmatch->outputmap(outImageName);
}
