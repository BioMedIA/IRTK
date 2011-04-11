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
	cerr << "-regionweight [value]              Region term weight for graphcut larger region term is more important"<<endl;
	cerr << "-timeweight [value]                Time weight for boundary term in 4D graphcut case larger time is more important" <<endl;
	cerr << "-outputname [name]                 Segmentation output prefix name" << endl;
	cerr << "-outputlabel [value]               Output Label in the output segmentation default 1" << endl;
	cerr << "-load                              Load segmentation results" << endl;
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
		connect = 0, outputvalue = 1,loadsegmentation = 0;
	double x,y,z,t,regionweight = 0.5,timeweight = 1.0;

	if( argc < 3 ) usage();

	// Number Of Images
	numberOfImages = atoi(argv[1]);
	argc--;
	argv++;
	irtkRealImage **input = new irtkRealImage *[numberOfImages];
	irtkRealImage **threshold = new irtkRealImage *[numberOfImages];
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
		if ((ok == false) && (strcmp(argv[1], "-load") == 0)) {
			argc--;
			argv++;
			loadsegmentation = 1;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-regionweight") == 0)) {
			argc--;
			argv++;
			regionweight = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-outputlabel") == 0)) {
			argc--;
			argv++;
			outputvalue = atoi(argv[1]);
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

		output[i] = new irtkRealImage;

		if(loadsegmentation)
			output[i]->Read(buffer);
		else{
			output[i]->Initialize(atr);
		}

		threshold[i] = new irtkRealImage;
		threshold[i]->Initialize(atr);

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
		cf.EvaluateGraphCut(threshold, input,atlas,numberOfImages,numberofatlas, timeweight,cutmode,connect,regionweight,0,c);
	}
	else{
		cerr << "multipleimage graphcut needs atlas information" << endl;
		exit(1);
	}
	for (n = 0; n < numberOfImages; n++) {
		for(int l = 0; l < threshold[n]->GetT(); l++){
			for(int k = 0; k < threshold[n]->GetZ(); k++){
				for(int j = 0; j < threshold[n]->GetY(); j++){
					for(int i = 0; i < threshold[n]->GetX(); i++){
						if(threshold[n]->GetAsDouble(i,j,k,l) > 0){
							output[n]->PutAsDouble(i,j,k,l,outputvalue);
						}
					}
				}
			}
		}
		if(!outputname)
			sprintf(buffer,"segtest%.2d.nii.gz",n);
		else
			sprintf(buffer,"%s%.2d.nii.gz",outputname,n);
		output[n]->Write(buffer);
	}

	for (i = 0; i < numberOfImages; i++) {
		delete input[i];
		delete threshold[i];
		delete output[i];
	}
	delete input;
	delete threshold;
	delete output;
	if(numberofatlas > 0){
		delete []atlas;
	}
}
