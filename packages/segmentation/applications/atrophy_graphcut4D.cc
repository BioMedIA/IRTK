#include <irtkGraphCutSegmentation_4D.h>
#include <irtkImage.h>
#include <string>
#include <irtkHistogram_1D.h>
#include <irtkDilation.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <irtkSubcorticalSegmentation_4D.h>

using namespace std;








void usage()
{
  cerr << "Usage: [N_timepoints] [image_1]...[image_N] [M_structures] [atlas_1]...[atlas_M] [wm_atlas] [gm_atlas] [csf_atlas] [lambda] [alpha] [output_dir]" << endl;
  cerr << "-minPerc [minimum percentage value for foreground model (default: 0.93)]" << endl;
  exit(1);
  
}

string getStructureName(string input){

	string str2 = "/";
	size_t found = 1000;
	int start = 0;
	while(found != string::npos){
		found = input.find(str2, start);
		if(found != string::npos)
			start = found+1;
	}
	if(start == 0){
		str2 = "\"";
	}
	while(found != string::npos){
		found = input.find(str2, start);
		if(found != string::npos)
			start = found+1;
	}
	size_t end = input.find(".", start);
	string output = input.substr(start,end-start);
	return output;
}

int main(int argc, char **argv){
	
	if (argc < 10) {
	    usage();
	    exit(1);
  	}



	int nrTimepoints = atoi(argv[1]);
	int nrStructures = atoi(argv[nrTimepoints+2]);
	string ** names = new string*[nrStructures];
	int ** structValues = new int*[nrStructures];

	cout << "Nr Timepoints: " << nrTimepoints << " Nr Structures: " << nrStructures << endl;

	if(argc < 8 + nrTimepoints + nrStructures) {
		usage();
		exit(1);
	}

	argc--;
	argv++;
	irtkRealImage tmp;
	cout << "Reading image " << argv[1] << endl;
	tmp.Read(argv[1]);
	string imageName = getStructureName(argv[1]);
	argc--;
	argv++;

	irtkImageAttributes attributes = tmp.GetImageAttributes();
	attributes._t = nrTimepoints;
	irtkRealImage * img = new irtkRealImage(attributes);

	for(int i = 0; i < nrTimepoints; i++){

		irtkRealPixel * ptr1 = tmp.GetPointerToVoxels(0, 0, 0);
		irtkRealPixel * ptr2 = img->GetPointerToVoxels(0, 0, 0, i);
		if(img->GetNumberOfVoxels()/nrTimepoints != tmp.GetNumberOfVoxels()){

			cerr << "Voxel dimensions do not agree!" << endl;
			exit(1);
		}
		else{
			for(int j = 0; j < tmp.GetNumberOfVoxels(); j++){
				*ptr2 = *ptr1;
				ptr1++;
				ptr2++;
			}
		}
		if(i != nrTimepoints-1){

			tmp.Read(argv[1]);
			cout << "Reading image " << argv[1] << endl;
			argc--;
			argv++;
		}
	}
	irtkRealImage ** structure_atlases;
	structure_atlases = new irtkRealImage*[nrStructures];
	argc--;
	argv++;
	for(int i = 0; i < nrStructures; i++){
		cout << "Reading atlas " << argv[1] << endl;
		string name (argv[1]);
		string s_name = getStructureName(name);
		names[i] = new string(s_name);
		structValues[i] = new int;
		*structValues[i] = i+1;

		structure_atlases[i] = new irtkRealImage(attributes);
		tmp.Read(argv[1]);
		for(int j = 0; j < nrTimepoints; j++){
			irtkRealPixel * ptr1 = tmp.GetPointerToVoxels(0, 0, 0);
			irtkRealPixel * ptr2 = structure_atlases[i]->GetPointerToVoxels(0, 0, 0, j);
			if(structure_atlases[i]->GetNumberOfVoxels()/nrTimepoints != tmp.GetNumberOfVoxels()){
				cerr << "Voxel dimensions do not agree!" << endl;
				exit(1);
			}
			else{
				for(int k = 0; k < tmp.GetNumberOfVoxels(); k++){
					*ptr2 = *ptr1;
					ptr1++;
					ptr2++;
				}

			}
		}
		argc--;
		argv++;
	}
	irtkRealImage ** tissue_priors = new irtkRealImage*[3];
	for(int i = 0; i < 3; i++){
		cout << "Reading tissue " << argv[1] << endl;
		tissue_priors[i] = new irtkRealImage(attributes);
		tmp.Read(argv[1]);
		for(int j = 0; j < nrTimepoints; j++){
			irtkRealPixel * ptr1 = tmp.GetPointerToVoxels(0, 0, 0);
			irtkRealPixel * ptr2 = tissue_priors[i]->GetPointerToVoxels(0, 0, 0, j);
			if(tissue_priors[i]->GetNumberOfVoxels()/nrTimepoints != tmp.GetNumberOfVoxels()){
				cerr << "Voxel dimensions do not agree!" << endl;
				exit(1);
			}
			else{
				for(int k = 0; k < tmp.GetNumberOfVoxels(); k++){
					*ptr2 = *ptr1;
					ptr1++;
					ptr2++;
				}

			}
		}
		argc--;
		argv++;
	}
	double lambda = atof(argv[1]);
	argc--;
	argv++;
	double alpha = atof(argv[1]);
	argc--;
	argv++;

	string outDir = argv[1];
	argc--;
	argv++;
	bool ok;
	double minPerc = 0.93;
	while (argc > 1){
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-minPerc") == 0)){
			argc--;
			argv++;
			minPerc = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}

		if (ok == false){
		  cerr << "Can not parse argument " << argv[1] << endl;
		  usage();
		}
	}
  	cout << "Min Percentage: " << minPerc << "Lambda: " << lambda << "Alpha: " << alpha << endl;
  	irtkSubcorticalSegmentation_4D * seg = new irtkSubcorticalSegmentation_4D();

  	cout<< endl << "init segmentation" << endl << endl;
  	seg->init(*img, structure_atlases, tissue_priors, nrStructures, outDir, names, structValues, imageName);
  	cout << "generate gaussians" << endl;
	seg->generateGaussians(minPerc);
	cout << "start segmentation" << endl;	

	double c = 0.5;
	double sigmaG = 0.02;

	seg->doSegmentation(lambda,alpha,c,sigmaG);
	

	

}






