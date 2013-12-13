/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/


#include <irtkImage.h>

#ifdef HAS_VTK

#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

// Default filenames
char *output_name = NULL, *prefix_name = NULL, *suffix_name = NULL, *textfile = NULL;

#define MAX_IMAGES 10000
#define EPSILON    0.001

void usage()
{
	cerr << "Usage: combinelandmarks [output] <input1..inputN> <options>" << endl;
	cerr << "" << endl;
	cerr << "Make an average landmark from a given set of inputs. All inputs must have the" << endl;
	cerr << "same number of landmarks. Names of input images can be given on the command line or" << endl;
	cerr << "in a file (see below)." << endl;
	cerr << "" << endl;
	cerr << "If they are named on the command line, landmarks are given equal weighting." << endl;
	cerr << "Different weights can be given if a file with names is provided." << endl;
	cerr << "" << endl;
	cerr << "where <options> is one or more of the following:" << endl;
	cerr << "<-names file>           File containing landmakr names and values to use for weighting." << endl;
	cerr << "                        The format should be one line for each file with" << endl;
	cerr << "                          \"name value\"" << endl;
	cerr << "                        on each line. If a mean and sigma value are specified" << endl;
	cerr << "<-prefix directory>     Directory in which to find the images." << endl;
	cerr << "<-scaling value>        Scaling value for final atlas (default = 1)\n";
	cerr << "<-maxnumber value>      Maxium number of images taken " << endl;
	cerr << "                          (default = off)\n";
	exit(1);
}

void finalize_landmarks(double *output, double scale, double *weight, int number)
{
	int i;
	double *ptr2,*ptr3;

	cout << "Constructing average landmarks...";
	cout.flush();
	ptr2 = output;
	ptr3 = weight;
	for (i = 0; i < number*3; i++) {
		if(*ptr3 > 0){
			*ptr2 = scale * *ptr2 / *ptr3;
		}else{
			*ptr2 = -1;
		}
		ptr2++;
		ptr3++;
	}
	cout << " done\n";
}

void add_landmark(vtkPoints *input, double *output, double inputweight, double *weight, int number)
{
	int i, j;
	double *ptr2, *ptr3;
	double p[3];

	cout << inputweight ;
	cout << " Adding landmark to combine..." << endl;
	cout.flush();

	if(!(input->GetNumberOfPoints() == number)){
		cerr << "Mismatch of number of landmarks:" << endl;
		exit(1);
	}

	ptr2 = output;
	ptr3 = weight;

	for (i = 0; i < number; i++) {

		input->GetPoint(i,p);

		for(j = 0; j < 3; j++){

			*ptr2 += p[j];
			*ptr3 += inputweight;
			ptr2++;
			ptr3++;
		}
	}

	cout << "done\n";
	cout.flush();
}

int main(int argc, char **argv)
{
	double *output = NULL;
	double *weight = NULL;
	double scale;
	int i, no, number, ok, maxiumimage = 10000;

	// Check command line
	if (argc < 4) {
		usage();
	}

	// Parse parameters
	output_name = argv[1];
	argc--;
	argv++;

	// Default: scaling factor
	scale = 1;

	// Parse input file names and values
	char **input_name = new char *[MAX_IMAGES];
	double *input_value = new double[MAX_IMAGES];
	double *input_weight = new double[MAX_IMAGES];
	vtkPolyDataReader *reader = NULL;
	reader = vtkPolyDataReader::New();
	vtkPolyData* surface = vtkPolyData::New();
	vtkPoints*   points  = NULL;

	// Parse any remaining paramters
	no = 0;

	while ((argc > 1) && (argv[1][0] != '-')) {
		input_name[no] = argv[1];
		input_weight[no] = 1;
		no++;
		argc--;
		argv++;
	}

	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-scaling") == 0)) {
			argc--;
			argv++;
			scale = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-names") == 0)) {
			argc--;
			argv++;
			textfile = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-maxnumber") == 0)) {
			argc--;
			argv++;
			maxiumimage = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-prefix") == 0)) {
			argc--;
			argv++;
			prefix_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-suffix") == 0)) {
			argc--;
			argv++;
			suffix_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if (ok == false) {
			cerr << "Unknown argument: " << argv[1] << endl << endl;
			usage();
		}
	}

	if (textfile != NULL) {
		ifstream in(textfile);
		if (in.is_open()) {
			while (!in.eof()) {
				input_name[no] = new char[256];
				in >> input_name[no] >> input_value[no];
				if (strlen(input_name[no]) > 0) {
					// Use input value directly as a weight.
					input_weight[no] = input_value[no];
					no++;
				}
			}
			in.close();
		} else {
			cout << "atlas: Unable to open file " << textfile << endl;
			exit(1);
		}
	}

	if(no > maxiumimage){
		no = maxiumimage;
	}

	// Read and add one input image after the other
	for (i = 0; i < no; i++) {
		char buffer[255];

		// Read input
		if (prefix_name != NULL) {
			sprintf(buffer, "%s%s", prefix_name, input_name[i]);
		} else {
			sprintf(buffer, "%s", input_name[i]);
		}
		input_name[i] = strdup(buffer);
		if (suffix_name != NULL) {
			sprintf(buffer, "%s%s", input_name[i], suffix_name);
		} else {
			sprintf(buffer, "%s", input_name[i]);
		}
		input_name[i] = strdup(buffer);

		cout << "Reading input landmark " << input_name[i] << "... ";
		cout.flush();
		reader->SetFileName(input_name[i]);
		reader->Update();
		
		surface->DeepCopy(reader->GetOutput());
		points = surface->GetPoints();
		
		if (i == 0) {
			number = points->GetNumberOfPoints();
			output = new double[number*3];
			for (int j = 0; j < number*3; j++)
				output[j] = 0;
			weight = new double[number*3];
			for (int j = 0; j < number*3; j++)
				weight[j] = 0;
		}

		add_landmark(points, output, input_weight[i], weight, number);

	}

	// Finalize atlas
	finalize_landmarks(output, scale, weight, number);

	// Write atlas
	cout << "Writing final landmarks to " << output_name << " ... " << endl;
	cout.flush();

	double p[3];
	
	double *ptr = output;
	for (int i=0; i < number; i++) {

		for (int j = 0; j < 3; j++){
			p[j] = *ptr;
			ptr++;
		}

		points->SetPoint(i,p);
	}

	// Write the final set
	vtkPolyDataWriter   *writer = vtkPolyDataWriter::New();
	writer->SetFileName(output_name);
	writer->SetInput(surface);
	writer->Write();

	cout << "done\n";
	cout.flush();
}

#else
#include <irtkCommon.h>
int main(int argc, char **argv)
{
	cerr << "combinelandmarks: this program needs to be compiled with vtk enabled.\n";
	return 0;
}
#endif