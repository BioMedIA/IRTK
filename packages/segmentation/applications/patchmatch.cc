/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkSegmentationFunction.h>

void usage()
{
	cerr << "Usage: patchmatch [target] [source] [patchradius] [outputmap]" << endl;
	cerr << "Where target and source are nD + t images, n indicates number of dimensionals of the image" << endl;
	cerr << "and t indicates the dimention of the normalized features at each pixel/voxel" << endl;
	cerr << "-searchradius         [0-1] random search radius in the image, 1 means whole image" << endl;
	cerr << "-nnfiterations        [N] number of iterations of the multi-atlas patchmatch" << endl;
	cerr << "-outputgraph filename output the distance graph between two image linked by the nnf" << endl;
	cerr << "-debug                open debug mode, output intermedia results" << endl;
	exit(1);

}
int main(int argc, char **argv){

	if (argc < 5) {
		usage();
	}

	char *output_name = NULL;
	char *output_graph_name = NULL;
	double xsize,ysize,zsize;
	int em_iterations = 1;
	int nnf_iterations = 40;
	double searchradius = 0.1;
	int debug = false;

	
	irtkGreyImage image;
	cout << "Reading target " << argv[1] << endl;
	image.Read(argv[1]);
	argc--;
	argv++;
	
	irtkGreyImage source;
	cout << "Reading source " << argv[1] << endl;
	source.Read(argv[1]);
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
		if ((ok == false) && (strcmp(argv[1], "-outputgraph") == 0)){
			argv++;
			argc--;
			output_graph_name = argv[1];
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

	irtkPatchMatch *patchmatch = new irtkPatchMatch(&image, &source, patchSize, 1);

	patchmatch->setDebug(debug);

	cout << "Creating patchmatch done"<<endl;

	cout << "Start optimization..."<<endl;

	patchmatch->setRandomrate(searchradius);

	patchmatch->run(nnf_iterations);

	cout << "Optimization done"<<endl;
		
	patchmatch->outputmap(outImageName);

	if(output_graph_name != NULL){
		patchmatch->outputgraph(output_graph_name);
	}
}