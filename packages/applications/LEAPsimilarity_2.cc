
#include <irtkRegistration.h>
#include <irtkPairwiseSimilarity.h>

#ifdef HAS_TBB
 #include "tbb/task_scheduler_init.h"
  #include "tbb/parallel_for.h"
 #include "tbb/blocked_range.h"


 using namespace tbb;


 class MultiThreadedLeapSimilarity{

 private:
	 int _nrImages;
	 int _sizeSet2;
	 bool _useMask;
	 irtkGreyImage ** _masks;
	 int _simType;
	 int _nrMasks;
	 bool _useSSD;
	 double * _results_nmi;
	// vector<pair<double,int> > * _results;
	 vector<string> _files;
	 string _imagedir;
	 irtkGreyImage * _target;




 public:
 	MultiThreadedLeapSimilarity(string imagedir, irtkGreyImage ** masks, double * results_nmi, irtkGreyImage *target, vector<string> files, bool useSSD){
 		_useMask = true; //always use a mask --> if not read in, it is generated here
		_simType = 1; //NMI
		_sizeSet2 = 1; // we only compare one image to all the others
		_nrMasks = 1; // always use one mask
		_nrImages = 1; // one image after the other

		_masks = masks;
		_results_nmi = results_nmi;
		_files = files;
		_imagedir = imagedir;
		_target = target;
		_useSSD = useSSD;
		if(_useSSD){
			_simType = 2; //SSD
		}


 	}

	void operator()(const blocked_range<int> &r) const {
		for (int i = r.begin(); i!= r.end(); i++) {


			irtkGreyImage ** _images;
			_images = new irtkGreyImage*[2]; // two images: current atlas and target image
			_images[0] = new irtkGreyImage(); // first image will be the current atlas
			_images[1] = new irtkGreyImage(); // second image will always be the target
			*_images[1] = *_target;


			cout << "Atlas " << i << ". Image ";

			irtkPairwiseSimilarity ps;
			ps.Initialize(_nrImages, _sizeSet2, _useMask, _masks, _simType, _nrMasks);
			string name = _imagedir + _files[i];
			_images[0]->Read(name.c_str());


			ps.LoadImages(_images);
			ps.GetSimilarities();
			double currentSim = ps.GetSimilarity(0,0,0);
			_results_nmi[i] = currentSim;
			delete _images[0];
			delete _images[1];
			delete [] _images;

		}
	}

 };

#endif


void usage()
{
	cerr << "usage LEAPsimilarity [atlasNames] [atlasDir] [targetImage] [ROI] [outFile]" << endl;
	cerr << "-SSD Use SSD instead of NMI as similarity measure." << endl;

		
}
int main(int argc, char **argv)
{
	if(argc < 6){
		usage();
		exit(1);
	}

	string imagefile = argv[1];
	argc--;
	argv++;
	string imagedir = argv[1];
	argc--;
	argv++;
	string targetImage = argv[1];
	argc--;
	argv++;
	string maskName = argv[1];
	argc--;
	argv++;
	string outname = argv[1];
	argc--;
	argv++;

	bool ok;
	bool useSSD = false;
	while (argc > 1){
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-SSD") == 0)){

		  argc--;
		  argv++;
		  useSSD = true;
		  ok = true;
		}
		if (ok == false){
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}
	//read imagenames
	vector<string> files = vector<string>();
	string line;
	ifstream from(imagefile.c_str());
	if (!from) {
		cerr << "Can't open file "
		<< endl;
		exit(1);
	}

	if (from.is_open())	{
		while (! from.eof() ) {
		  getline (from,line);
		  if(line != ""){
		  	string fileName = imagedir+ line;
		  	files.push_back(line);
		  }
		}
	from.close();
	}

	//nrImages gives the number of atlas images excluding the target.
//	int nrImages = files.size();
//	files.push_back(targetImage);

//	cout << "Nr of regions: " << endl;


	irtkGreyImage ** masks;
	masks = new irtkGreyImage*[1];
	masks[0] = new irtkGreyImage();
	masks[0]->Read(maskName.c_str());
	irtkGreyPixel *ptr = masks[0]->GetPointerToVoxels();
	for(int i = 0; i < masks[0]->GetNumberOfVoxels(); i++){
		*ptr = *ptr - 1;
		ptr++;
	}

	double *results = new double[files.size()];
	irtkGreyImage * target;
	target = new irtkGreyImage();
	target->Read(targetImage.c_str());

	#ifdef HAS_TBB
	task_scheduler_init init;
	MultiThreadedLeapSimilarity evaluate(imagedir, masks, results, target, files, useSSD);
	int blocks = 2;
	parallel_for(blocked_range<int>(0, int(files.size()), int(blocks)), evaluate);

	#else

	bool useMask = true; //always use a mask
	int simType = 1; //NMI
	out << "Using NMI as similarity metric" << endl;
	if(useSSD){
		simType = 2; //SSD
		cout << "Using SSD as similarity metric" << endl;
	}
	int sizeSet2 = 1; // we only compare one image to all the others
	int nrMasks = 1; // always use one mask
	int nrImages = 1; // one image after the other



	irtkGreyImage ** images;
	images = new irtkGreyImage*[2]; // two images: current atlas and target image
	images[0] = new irtkGreyImage(); // first image will be the current atlas
	images[1] = new irtkGreyImage(); // second image will always be the target
	for(int i = 0; i < int(files.size()); i++){


		       *images[1] = *target;

		cout << "Atlas " << i << ". Image ";

		irtkPairwiseSimilarity ps;

		ps.Initialize(nrImages, sizeSet2, useMask, masks, simType, nrMasks);
		string name = imagedir + files[i];
		images[0]->Read(name.c_str());
		ps.LoadImages(images);
		ps.GetSimilarities();
		double currentSim = ps.GetSimilarity(0,0,0);
		results[i] = currentSim;
	}

	#endif

	cout << "Writing " << outname << endl;
	ofstream output;
	ostringstream stm;
	output.open(outname.c_str());
	for(int i = 0; i < (int)files.size(); i++){
		output << results[i] << endl;
	}
	output.close();


    delete []results;

}
