

#include <irtkPatchBasedSegmentation.h>
#include <irtkNormalizeNyul.h>





using namespace std;

void usage()
{
  cerr << "Usage: [N] [imageFile] [atlasDir] [labelDir] [dofDir] [imageName] [atlasID1,...,atlasIDN] [patchSize] [neighbourhoodsize] [output]" << endl;
  cerr << "-mask" << endl;
  cerr << "-normMask" << endl;
  cerr << "-winnerTakesAll" << endl;
  exit(1);
  
}
int main(int argc, char **argv){

	if (argc < 8) {
		    usage();
		    exit(1);
	  }
	int verbose = 1;
	int nAtlases = atoi(argv[1]);
	argv++;
	argc--;
	int padding = -1;
    
	if (argc < nAtlases + 8) {
	    usage();
	    exit(1);
  	}

	string imagefilename = argv[1];
	argc--;
	argv++;

	vector<string> files = vector<string>();
	string line;
	ifstream from(imagefilename.c_str());
	if (!from) {
		cerr << "Can't open file "
		<< endl;
		exit(1);
	}

	if (from.is_open())	{
		while (! from.eof() ) {
		  getline (from,line);
		  string fileName = line;
		  files.push_back(fileName);
		}
	from.close();
	}

	string atlasDir = argv[1];
	argv++;
	argc--;
	string segDir = argv[1];
	argv++;
	argc--;
	string dofDir = argv[1];
	argv++;
	argc--;
	string imageName = argv[1];
	argv++;
	argc--;

	irtkRealImage image;
	string imageDirName = atlasDir + imageName;
	image.Read(imageDirName.c_str());



	irtkRealImage ** atlases;
	irtkGreyImage ** labels;

	atlases = new irtkRealImage*[nAtlases];
	labels = new irtkGreyImage*[nAtlases];
	int ids[nAtlases];
	for(int i = 0; i < nAtlases; i++){
		ids[i] = atoi(argv[1]);
		argv++;
		argc--;
	}

	int patchSize = atoi(argv[1]);
	argc--;
	argv++;
	patchSize = patchSize / 2;
	cout << "Patch size: " << patchSize << endl;
	patchSize = floor(patchSize);

	int neighbourhoodsize = atoi(argv[1]);
	argc--;
	argv++;
	neighbourhoodsize = neighbourhoodsize / 2;
	cout << "Neighbourhood size: " << neighbourhoodsize << endl;
	neighbourhoodsize = floor(neighbourhoodsize);

	char * outName = argv[1];
	argc--;
	argv++;
	bool ok;
	irtkGreyImage mask, normMask;
	bool useMask = false;
	bool useNormMask = false;
	bool useManifold = false;
	bool winnerTakesAll = false;
	while (argc > 1){
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
			argc--;
			argv++;
			useMask = true;
			mask.Read(argv[1]);
			argc--;
			argv++;

			ok = true;
		}

		if ((ok == false) && (strcmp(argv[1], "-normMask") == 0)){
			argc--;
			argv++;
			useNormMask = true;
			normMask.Read(argv[1]);
			argc--;
			argv++;

			ok = true;
		}

		if ((ok == false) && (strcmp(argv[1], "-winnerTakesAll") == 0)){
			argc--;
			argv++;
			winnerTakesAll = true;

			ok = true;
		}

		if (ok == false){
		  cerr << "Can not parse argument " << argv[1] << endl;
		  usage();
		}
	}


	// Finished reading parameters

	if(useNormMask){
		irtkRealPixel * ptr = image.GetPointerToVoxels();
		irtkGreyPixel * mPtr = normMask.GetPointerToVoxels();
		for(int i = 0; i < image.GetNumberOfVoxels(); i++){
			if(*mPtr<1){
				*ptr = 0;
			}
			ptr++;
			mPtr++;
		}
	}

	image.PutMinMax(0, 100);

	if(useNormMask){
		irtkRealPixel * ptr = image.GetPointerToVoxels();
		irtkGreyPixel * mPtr = normMask.GetPointerToVoxels();
		for(int i = 0; i < image.GetNumberOfVoxels(); i++){
			if(*mPtr<1){
				*ptr = padding;
			}
			ptr++;
			mPtr++;
		}
	}

	if(verbose)
		cout << "Reading atlases, transforming to target space, normalizing intensities" << endl;


	for(int i = 0; i < nAtlases; i++){

		if(verbose)
			cout << files[ids[i]] << endl;

		irtkGreyImage atlasInput;
		string filename = atlasDir + files[ids[i]];
		atlasInput.Read(filename.c_str());

		atlases[i] = new irtkRealImage;

		atlases[i]->Initialize(image.GetImageAttributes());

		irtkGreyImage labelInput;
		filename = segDir + files[ids[i]];
		labelInput.Read(filename.c_str());


		labels[i] = new irtkGreyImage;
		labels[i]->Initialize(image.GetImageAttributes());

		/////
		string atlasName = files[ids[i]];
		int pos = atlasName.find('.nii.gz');
		string atlasBaseName = atlasName.substr(0,pos-6);
		string dofNameAtlas =  dofDir + "/areg-" + atlasBaseName + ".dof.gz" ;
		char* pch = (char*)malloc( sizeof( char ) *(dofNameAtlas.length() +1) );
		strcpy( pch, dofNameAtlas.c_str() );

		irtkTransformation *transformAtlas = irtkTransformation::New(pch);
		irtkAffineTransformation *affineAtlas = dynamic_cast<irtkAffineTransformation *>(transformAtlas);

		pos = imageName.find('.nii.gz');
		string imageBaseName = imageName.substr(0,pos-6);
		string dofNameImage =  dofDir + "/areg-" + imageBaseName + ".dof.gz" ;
		pch = (char*)malloc( sizeof( char ) *(dofNameImage.length() +1) );
		strcpy( pch, dofNameImage.c_str() );


		irtkTransformation *transformImage = irtkTransformation::New(pch);
		irtkAffineTransformation *affineImage = dynamic_cast<irtkAffineTransformation *>(transformImage);
		affineImage->Invert();

		// Convert to matrices
		irtkMatrix matrix1 = affineAtlas->GetMatrix();
		irtkMatrix matrix2 = affineImage->GetMatrix();
		irtkMatrix matrixOut = matrix1 * matrix2;
		irtkAffineTransformation  affineComposedTrans;
		affineComposedTrans.PutMatrix(matrixOut);
		irtkImageFunction *interpolator = new irtkLinearInterpolateImageFunction;
		irtkImageFunction *interpolatorNN = new irtkNearestNeighborInterpolateImageFunction;
		irtkImageTransformation *imagetransformation = new irtkImageTransformation;
		imagetransformation->SetInput(&atlasInput, &affineComposedTrans);
		imagetransformation->SetOutput(atlases[i]);
		imagetransformation->PutInterpolator(interpolator);
		cout << "Transforming atlas " << "...";

		imagetransformation->Run();
		cout << "done" << endl;
		imagetransformation->SetInput(&labelInput, &affineComposedTrans);
		imagetransformation->SetOutput(labels[i]);
		imagetransformation->PutInterpolator(interpolatorNN);
		cout << "Transforming label " << "...";
		imagetransformation->Run();
		cout << "done" << endl;
		delete interpolator;
		delete imagetransformation;

		//pad atlas

		if(useNormMask){
			irtkRealPixel * ptr = atlases[i]->GetPointerToVoxels();
			irtkGreyPixel * mPtr = normMask.GetPointerToVoxels();
			for(int j = 0; j < atlases[i]->GetNumberOfVoxels(); j++){
				if(*mPtr<1){
					*ptr = 0;
				}
				ptr++;
				mPtr++;
			}
		}
		//Do rescaling 0-100

		atlases[i]->PutMinMax(0, 100);

		if(useNormMask){
			irtkRealPixel * ptr = atlases[i]->GetPointerToVoxels();
			irtkGreyPixel * mPtr = normMask.GetPointerToVoxels();
			for(int j = 0; j < atlases[i]->GetNumberOfVoxels(); j++){
				if(*mPtr<1){
					*ptr = padding;
				}
				ptr++;
				mPtr++;
			}
		}

		/////Do normalization

		cout << "Start normalization...";
		irtkNormalizeNyul nn(*atlases[i], image);
		nn.SetPadding(padding, padding);
		nn.Run();
		*atlases[i] = nn.GetOutput();
		cout << "done" << endl;



	}



	irtkPatchBasedSegmentation patchBased(image, atlases, labels, nAtlases, patchSize, neighbourhoodsize);
	if(useMask)
		patchBased.SetMask(mask);
	if(winnerTakesAll){
		cout << "Winner takes all" << endl;
		patchBased.WinnerTakesAll();
	}
	patchBased.SetPadding(padding);
	cout << "Run pb" << endl;
	patchBased.Run();
	patchBased.GetConsensusSegmentation();
	cout << "Writing output" << endl;
	patchBased.WriteSegmentation(outName);



}



