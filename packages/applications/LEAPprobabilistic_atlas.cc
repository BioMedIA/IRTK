
#include <irtkRegistration.h>

irtkMatrix read_csv(string csvFilename1, string csvFilename2);
void writeAtlases(string filename, int * atlases, int nrAtlases);

void usage()
{
	cerr << "usage LEAPPropabilistic_atlas [N] [segmentation1,...,segmentationN] [label] [output]" << endl;

		
}
int main(int argc, char **argv)
{

	if(argc < 4){
		usage();
		exit(1);
	}

	int N = atoi(argv[1]);
	argc--;
	argv++;
	if(argc != N+3){
		usage();
		exit(1);
	}
	irtkGreyImage ** input;
	input = new irtkGreyImage*[N];
	for(int i = 0; i < N; i++){
		input[i] = new irtkGreyImage;
		input[i]->Read(argv[1]);
		argc--;
		argv++;
	}
	int label = atoi(argv[1]);
	argc--;
	argv++;
	irtkRealImage output;
	output.Initialize(input[0]->GetImageAttributes());
	string outname = argv[1];
	argc--;
	argv++;
	double max = 0;
	for(int i = 0; i < N; i++){
		irtkRealPixel *ptr1 = output.GetPointerToVoxels();
		irtkGreyPixel *ptr2 = input[i]->GetPointerToVoxels();
		for(int j = 0; j < input[0]->GetNumberOfVoxels(); j++){
			if(*ptr2==label)
				*ptr1 += double(1);
			if(*ptr1 > max)
				max = *ptr1;
			ptr1++;
			ptr2++;
		}

	}
	cout << "Maxval: " << max << endl;
	irtkRealPixel *ptr1 = output.GetPointerToVoxels();
	for(int i = 0; i < input[0]->GetNumberOfVoxels(); i++){


		if(*ptr1>0){
			*ptr1 /= max;
			*ptr1 *= 100;
		}
		ptr1++;
	}

	output.Write(outname.c_str());

}







