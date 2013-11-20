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
	cerr << "Usage: mapmSegmentation [N] [atlasfile] [target] [labelfile] [outputlabel] [patchradius in mm]" << endl;
	cerr << "-labeldistances       [targetlabeldistanceimage] [atlaslabeldistancefile]" << endl;
	cerr << "-searchradius         [0-1] random search radius in the image, 1 means whole image" << endl;
	cerr << "-nnfiterations        [N] number of iterations of the multi-atlas patchmatch" << endl;
	cerr << "-emiterations         [N] number of iterations of EM algorithm" << endl;
	cerr << "-output               [filename] output final mapping to the file" << endl;
	cerr << "-distanceweight       [weight] weight for the label distances" << endl;
	cerr << "-directionalradius    [x y z] radius for each direction in mm" << endl;
	cerr << "-debug                open debug mode, output intermedia results" << endl;
	cerr << "where atlasfile, labelfile and atlaslabeldistancefile are the txt file containing the location of the images" << endl;
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
	irtkGreyImage ** labels = NULL;
	irtkRealImage * targetdistance = NULL;
	irtkRealImage ** atlasdistances = NULL;
	char *output_name = NULL;
	int em_iterations = 1;
	int nnf_iterations = 40;
	double searchradius = 0.1;
	double weight = 1;
	double radius_x = 0;
	double radius_y = 0;
	double radius_z = 0;
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

	char * outLabelName = argv[1];
	argc--;
	argv++;

	int patchSize = atoi(argv[1]);
	argc--;
	argv++;
	cout << "Patch radius: " << patchSize << endl;

	while (argc > 1){
		int ok = false;
		if ((ok == false) && (strcmp(argv[1], "-labeldistances") == 0)){
			argv++;
			argc--;
			targetdistance = new irtkRealImage(argv[1]);
			argv++;
			argc--;
			atlasdistances = new irtkRealImage*[nAtlases];

			line = "";
			infile.open (argv[1]);

			argv++;
			argc--;

			for(int i = 0; i < nAtlases; i++){
				if(!infile.eof()){
					getline(infile,line);
					cout << "Reading label distance images" << line << endl;
					atlasdistances[i] = new irtkRealImage((char*)line.c_str());
				}else{
					cout << "Not enough label distance images, should be " << nAtlases << " actually " << i+1 << endl;
					break;
				}
			}

			infile.close();
			ok = true;
		}
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
		if ((ok == false) && (strcmp(argv[1], "-directionalradius") == 0)){
			argv++;
			argc--;
			radius_x = atof(argv[1]);
			argv++;
			argc--;
			radius_y = atof(argv[1]);
			argv++;
			argc--;
			radius_z = atof(argv[1]);
			argv++;
			argc--;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-distanceweight") == 0)){
			argv++;
			argc--;
			weight = atof(argv[1]);
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

	irtkGreyImage label;

	irtkMAPatchMatchSegmentation * patchmatchesegmentation = 
		new irtkMAPatchMatchSegmentation(&image,atlases,targetdistance,atlasdistances, &label, labels, patchSize,nAtlases);

	patchmatchesegmentation->setDebug(debug);

	cout << "Creating patchmatch done"<<endl;

	cout << "Start optimization..."<<endl;


	patchmatchesegmentation->setRandomrate(searchradius);

	patchmatchesegmentation->setWeight(weight);

	if(radius_x > 0 && radius_y > 0 && radius_z > 0){
		patchmatchesegmentation->setDirectionalRadius(radius_x, radius_y, radius_z);
	}

	patchmatchesegmentation->run(nnf_iterations);

	cout << "Optimization done"<<endl;

	patchmatchesegmentation->generateLabels();

	cout << "Writing output" << endl;
	label.Write(outLabelName);
	
	if(output_name != NULL){
		patchmatchesegmentation->outputmap(output_name);
	}

}



