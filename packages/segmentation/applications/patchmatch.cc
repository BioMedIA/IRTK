/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkPatchMatch.h>

void usage()
{
	cerr << "Usage: patchmatch [N] [atlasfile] [image] [patchradius] [outputimage]" << endl;
	cerr << "-labels               [labelfile] [outputlabel]" << endl;
	cerr << "-levels               [N] use multiple level similarity between patches" << endl;
	cerr << "-LDimage              [image] reference low resolution image" << endl;
	cerr << "-searchradius         [0-1] random search radius in the image, 1 means whole image" << endl;
	cerr << "-nnfiterations        [N] number of iterations of the multi-atlas patchmatch" << endl;
	cerr << "-emiterations         [N] number of iterations of EM algorithm" << endl;
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
	irtkGreyImage * ld_image = NULL;
	irtkGreyImage ** labels = NULL;
	double xsize,ysize,zsize;
	int levels;
	int slevels = 1;
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
		if ((ok == false) && (strcmp(argv[1], "-labels") == 0)){
			argv++;
			argc--;
			labels = new irtkGreyImage*[nAtlases];

			line = "";
			infile.open (argv[1]);

			argv++;
			argc--;

			for(int i = 0; i < nAtlases; i++){
				if(!infile.eof()){
					getline(infile,line);
					cout << "Reading labels " << line << endl;
					labels[i] = new irtkGreyImage((char*)line.c_str());
				}else{
					cout << "Not enough labels, should be " << nAtlases << " actually " << i+1 << endl;
					break;
				}
			}

			infile.close();

			outLabelName = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-LDimage") == 0)){
			argv++;
			argc--;
			ld_image = new irtkGreyImage(argv[1]);
			argv++;
			argc--;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-levels") == 0)){
			argv++;
			argc--;
			slevels = atoi(argv[1]);
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
		if ((ok == false) && (strcmp(argv[1], "-emiterations") == 0)){
			argv++;
			argc--;
			em_iterations = atoi(argv[1]);
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

	irtkImageFunction *interpolator = new irtkBSplineInterpolateImageFunction;
	xsize = atlases[0]->GetX();
	ysize = atlases[0]->GetY();
	zsize = atlases[0]->GetZ();

	levels = 0;
	/*while(xsize > patchSize*4.0 
		&& ysize > patchSize*4.0 
		&& zsize > patchSize*4.0){
		levels++;
		xsize /= 2;
		ysize /= 2;
		zsize /= 2;
	}*/

	if(levels < 1)
		levels = 1;

	cout << "Number of levels: " << levels << endl;

	atlases[0]->GetPixelSize(&xsize, &ysize, &zsize);

	irtkPatchMatch ** patchmatches = new irtkPatchMatch*[levels];
	irtkGreyImage *** latlases = new irtkGreyImage**[levels];
	irtkGreyImage ** limages = new irtkGreyImage*[levels];

	cout << "Creating pymaids..."<<endl;

	for(int i = 0; i < levels; i++){

		cout << "level: " << i << endl;

		limages[i] = new irtkGreyImage(image);

		irtkResampling<short> resampling(xsize, ysize, zsize);
		resampling.SetInput (limages[i]);
		resampling.SetOutput(limages[i]);
		resampling.SetInterpolator(interpolator);
		resampling.Run();

		latlases[i] = new irtkGreyImage*[nAtlases];

		for(int j = 0; j < nAtlases; j++){

			latlases[i][j] = new irtkGreyImage(*atlases[j]);

			irtkResampling<short> resampling(xsize, ysize, zsize);
			resampling.SetInput (latlases[i][j]);
			resampling.SetOutput(latlases[i][j]);
			resampling.SetInterpolator(interpolator);
			resampling.Run();   
		}

		patchmatches[i] = new irtkPatchMatch(limages[i], latlases[i], patchSize, nAtlases, slevels);

		if(i == 0 && ld_image != NULL){
			patchmatches[i]->setDecimatedImage(ld_image);
			em_iterations++;
		}

		patchmatches[i]->setDebug(debug);

		xsize *= 2;
		ysize *= 2;
		zsize *= 2;
	}

	cout << "Creating pymaids done"<<endl;

	cout << "Start optimization..."<<endl;

	for(int i = levels - 1; i >= 0 ; i--){
		cout << "level: " << i << endl;

		patchmatches[i]->setRandomrate(searchradius);

		for(int j = 0; j < em_iterations; j++){
		
		  patchmatches[i]->runEMIteration(nnf_iterations);
		
		}

		if(i > 0){
			cout << "upsample fields to level: " << i-1 << endl;
			patchmatches[i]->upSampleNNF(patchmatches[i-1]);
		}
	}

	cout << "Optimization done"<<endl;

	cout << "Writing output" << endl;
	limages[0]->Write(outImageName);

	///write label segmentation;
	if(outLabelName != NULL){
		irtkGreyImage label;
		patchmatches[0]->generateLabels(&label, labels);
		label.Write(outLabelName);
	}

}



