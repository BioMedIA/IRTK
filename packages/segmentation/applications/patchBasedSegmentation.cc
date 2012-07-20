#include <irtkPatchBasedSegmentation.h>

using namespace std;

void usage()
{
  cerr << "Usage: [N] [atlas1,...,atlasN] [labels1,...labelsN] [image] [patchSize] [neighbourhoodsize] [output]" << endl;
  cerr << "-mask" << endl;
  cerr << "-probabilistic" << endl;
  cerr << "-padding" << endl;
  exit(1);
  
}
int main(int argc, char **argv){

	if (argc < 6) {
		    usage();
		    exit(1);
	  }

	int nAtlases = atoi(argv[1]);
	cout << "N Atlases: " << nAtlases << endl;
	argv++;
	argc--;
    
	if (argc < nAtlases*2 + 4) {
	    usage();
	    exit(1);
  	}

	irtkRealImage ** atlases;
	irtkGreyImage ** labels;

	atlases = new irtkRealImage*[nAtlases];
	labels = new irtkGreyImage*[nAtlases];

	for(int i = 0; i < nAtlases; i++){
		atlases[i] = new irtkRealImage;
		atlases[i]->Read(argv[1]);
		argv++;
		argc--;
	}

	for(int i = 0; i < nAtlases; i++){
		labels[i] = new irtkGreyImage;
		labels[i]->Read(argv[1]);
		argv++;
		argc--;
	}

	irtkRealImage image;
	cout << "Reading image " << argv[1] << endl;
	image.Read(argv[1]);
	argc--;
	argv++;
	int patchSize = atoi(argv[1]);
	argc--;
	argv++;
	patchSize = patchSize / 2;
	cout << "Patch size: " << patchSize << endl;


	int neighbourhoodsize = atoi(argv[1]);
	argc--;
	argv++;
	neighbourhoodsize = neighbourhoodsize / 2;
	cout << "Neighbourhood size: " << neighbourhoodsize << endl;

	char * outName = argv[1];
	argc--;
	argv++;
	bool ok;
	irtkGreyImage mask;
	bool useMask = false;
	bool useProb = false;
	int padding = -1;
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

		if ((ok == false) && (strcmp(argv[1], "-probabilistic") == 0)){
			argc--;
			argv++;
			useProb = true;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-padding") == 0)){
			argc--;
			argv++;
			padding = atoi(argv[1]);
			ok = true;
		}


		if (ok == false){
		  cerr << "Can not parse argument " << argv[1] << endl;
		  usage();
		}
	}

	irtkPatchBasedSegmentation patchBased(image, atlases, labels, nAtlases, patchSize, neighbourhoodsize);
	if(useMask)
		patchBased.SetMask(mask);
	patchBased.SetPadding(padding);
	patchBased.Run();
	if(useProb)
		patchBased.GetProbabilisticSegmentation();
	else
		patchBased.GetConsensusSegmentation();
	cout << "Writing output" << endl;
	patchBased.WriteSegmentation(outName);



}



