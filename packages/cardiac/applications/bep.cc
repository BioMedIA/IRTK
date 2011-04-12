#include <irtkBep.h>

int Zcheck = 0;

void usage()
{
	cerr << "Usage: bep [originalsequence] [transformation] [PointSet] [ResultOutput] [MaxOutput] [SegmentationOutput]" << endl;
	cerr << "where PointSet is the sector landmarks of the bull's eye plot" << endl;
	cerr << "AA is the slice landmark between apex and apical" << endl;
	cerr << "AM is the slice landmark between apical and mid" << endl;
	cerr << "MB is the slice landmark between mid and basal" << endl;
	cerr << "-threshold [thresholdname] segmentation information for segmentation"    << endl;
	cerr << "-segmentation [thresholdname] segmentation information for segmentation"    << endl;
	cerr << "-dofprefix [value] transformation filename prefix"    << endl;
	cerr << "-mode [0/1/2/3/4/5] radial/motion/strain/longitudinal strain/ radial and circumferential strain/normalized strain" << endl;
	exit(1);
}

void checkslice(irtkGreyImage& target, int z){
	if(z < 0 || z > target.GetZ()){
		cerr<<"Invaliade z slice number: "<< z <<endl;
		exit(1);
	}
	if(z == target.GetZ() && Zcheck == 1){
		cerr<<"Invaliade z slice number: "<< z <<endl;
		exit(1);
	}
	return;
}

int main( int argc, char** argv )
{
	irtkGreyImage target,source,late;
	irtkPointSet landmarks;
	irtkPointSet olandmarks;
	int tip,aa,am,mb,bottom,swapped,mode,prefix;
	double gap = 0;
	char *outputfilename = NULL;
	char *maxfilename = NULL;
	char *osequence = NULL;
	char *segout = NULL;
	char *thresholdname = NULL;
	char *segmentname = NULL;
	irtkRealImage threshold,compare,segmentation;
	irtkImageAttributes atr;
	irtkBep bf;
	int i,j,k,l,ok;

	if( argc < 7 ) usage();
	
	tip = 0;
	mode = 0;
	prefix = -1;
	aa = 0; am = 0; mb = 0; bottom = 0;

	// Read image
	cout << "Reading original image: " << argv[1] << endl;
	target.Read(argv[1]);
	argc--;
	argv++;
	osequence = argv[1];
	argc--;
	argv++;
	olandmarks.ReadVTK(argv[1]);
	argc--;
	argv++;
	outputfilename = argv[1];
	argc--;
	argv++;
	maxfilename = argv[1];
	argc--;
	argv++;
	segout = argv[1];
	argc--;
	argv++;
	remove(outputfilename);
	remove(maxfilename);

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-threshold") == 0)) {
			argc--;
			argv++;
			thresholdname = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-segmentation") == 0)) {
			argc--;
			argv++;
			segmentname = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-mode") == 0)) {
			argc--;
			argv++;
			mode = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-prefix") == 0)) {
			argc--;
			argv++;
			prefix = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	if(!thresholdname){
		//find out number of time frame
		cerr << "Segmentation result not given please use graphcut4D to generate segmentation"<< endl;
		exit(1);
	}else{
		threshold.Read(thresholdname);
	}

	//Initialize landmark
	if(!segmentname){
		for( i=0; i<olandmarks.Size(); i++){
			target.WorldToImage(olandmarks(i));
		}
		aa = round(olandmarks(5)._z);
		bottom = round(olandmarks(6)._z);
		swapped = 0;
		if(aa>bottom){
			swap(aa,bottom);
			swapped = 1;
		}
		l = 0;
		for( k = aa; k < bottom; k++){
			for( j = 0; j < threshold.GetY(); j++){
				for( i = 0; i < threshold.GetX(); i++){
					if(threshold.GetAsDouble(i,j,k) == 2 && l == 0){
						l = k;
					}
				}
			}
		}
		aa = l;
		gap = (double)(bottom - aa)/3.0;
		am = round(aa + gap);
		mb = round(bottom - gap); 	
		checkslice(target,aa);
		checkslice(target,am);
		checkslice(target,mb);
		checkslice(target,bottom);
		if(swapped == 1){
			swap(aa,bottom);
			swap(am,mb);
		}
		tip = 0;
	}

	bf.SetInput(target,threshold,osequence);
	if(!segmentname){
		bf.SetLandmarks(olandmarks,tip,aa,am,mb,bottom);
	}
	bf.SetOutput(outputfilename,maxfilename);

	//Generate PointSet from 3 landmarks
	if(!segmentname){
		bf.Initialize();
		bf.GenerateSegmentation(segout);
	}else{
		segmentation.Read(segmentname);
		bf.SetSegmentation(segmentation);
	}

	//Now start
	if(mode == 0)
		bf.AnalysisRadius(prefix);
	else if(mode == 1){
		if(prefix == -1)
			bf.AnalysisMotion();
		else
			bf.AnalysisMotion(prefix);
	}else if(mode >= 2){
		if(prefix == -1)
			bf.AnalysisStrain(3,mode - 2);
		else
			bf.AnalysisStrain(prefix, mode - 2);
	}
	bf.Bullseyeplot();
}