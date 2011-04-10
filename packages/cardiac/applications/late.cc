#include <irtkBep.h>

int Zcheck = 0;

void usage()
{
	cerr << "Usage: late [originalsequence] [lateimage] [PointSet] [Output] [MaxOutput]" << endl;
	cerr << "-threshold [thresholdname] segmentation information for segmentation"    << endl;
	cerr << "-mod toggle late/wall/motion -1/0/1" << endl;
	cerr << "where PointSet is the sector landmarks of the bull's eye plot" << endl;
	exit(1);
}

void checkslice(irtkGreyImage& target, int z){
	if(z < 0 || z > target.GetZ()){
		cerr<<"Invaliade z slice number: "<< z <<endl;
		exit(1);
	}
	if(z == target.GetZ() && Zcheck == 0){
		Zcheck = 1;
	}
	if(z == target.GetZ() && Zcheck == 1){
		cerr<<"Invaliade z slice number: "<< z <<endl;
		exit(1);
	}
	return;
}

int main( int argc, char** argv )
{
	irtkGreyImage target,late;
	irtkPointSet landmarks;
	irtkPointSet olandmarks;
	int tip,aa,am,mb,bottom;
	double gap;
	char *outputfilename1 = NULL;
	char *outputfilename2 = NULL;
	char *osequence = NULL;
	char *lsequence = NULL;
	char *thresholdname = NULL;
	char *segmentname = NULL;
	irtkRealImage threshold,compare,segmentation;
	irtkImageAttributes atr;
	irtkBep bf;
	int i,j,k,l,t,ok, latewall = 1,swapped;

	if( argc < 6 ) usage();

	// Read image
	cout << "Reading image: " << argv[1] << endl;
	target.Read(argv[1]);
	argc--;
	argv++;
	cout << "Reading late image: " << argv[1] << endl;
	compare.Read(argv[1]);
	argc--;
	argv++;
	olandmarks.ReadVTK(argv[1]);
	argc--;
	argv++;
	outputfilename1 = argv[1];
	argc--;
	argv++;
	outputfilename2 = argv[1];
	argc--;
	argv++;
	remove(outputfilename1);
	remove(outputfilename2);

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-mode") == 0)) {
			argc--;
			argv++;
			latewall = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
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
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}  

	//find out number of time frame
	atr = target.GetImageAttributes();
	t = atr._t;
	atr._t = 1;
	if(!thresholdname){
		cerr << "Segmentation result not given please use graphcut4D to generate segmentation"<< endl;
	}else{
		threshold.Read(thresholdname);
	}

	//Initialize landmark
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

	bf.SetInput(target,target,threshold,osequence,lsequence);
	bf.SetLandmarks(olandmarks,tip,aa,am,mb,bottom);
	bf.SetOutput(outputfilename1,outputfilename2);
	if(latewall == 0){
		for(k = 0; k < atr._z; k++){
			for(j = 0; j < atr._y; j++){
				for(i = 0; i < atr._x; i++){
					if(compare.GetAsDouble(i,j,k) > 0){
						compare.PutAsDouble(i,j,k,1);
					}
				}
			}
		}	
		bf.SetCompare(compare);		
	}else{
		bf.SetCompare(compare);
	}

	//Generate PointSet from 3 landmarks
	if(!segmentname){
	    bf.Initialize();
		bf.GenerateSegmentation("surface");
	}else{
		segmentation.Read(segmentname);
		bf.SetSegmentation(segmentation);
	}
	//start
	bf.Bullseyeplot(latewall);

	bf.Finalize();
}
