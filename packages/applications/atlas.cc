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

#include <irtkFileToImage.h>

#include <irtkHistogram.h>

#include <irtkTransformation.h>

// Default filenames
char *output_name = NULL, *prefix_name = NULL, *suffix_name = NULL, *textfile = NULL;

#define MAX_IMAGES 10000
#define EPSILON    0.001

void usage()
{
	cerr << "Usage: atlas [output] <input1..inputN> <options>" << endl;
	cerr << "" << endl;
	cerr << "Make an atlas from a given set of input images. All input images must have the" << endl;
	cerr << "same voxel lattice. Names of input images can be given on the command line or" << endl;
	cerr << "in a file (see below)." << endl;
	cerr << "" << endl;
	cerr << "If they are named on the command line, images are given equal weighting." << endl;
	cerr << "Different weights can be given if a file with names is provided." << endl;
	cerr << "" << endl;
	cerr << "where <options> is one or more of the following:" << endl;
	cerr << "<-imagenames file>      File containing image names and values to use for weighting." << endl;
	cerr << "                        The format should be one line for each image with" << endl;
	cerr << "                          \"image_name value\"" << endl;
	cerr << "                        on each line. If a mean and sigma value are specified" << endl;
	cerr << "                        with the \'-gaussian\' flag (see below), then the values" << endl;
	cerr << "                        are converted to weights using the corresponding Gaussian" << endl;
	cerr << "                        function. Otherwise, the values in the file are used" << endl;
	cerr << "                        directly as weights." << endl;
	cerr << "<-prefix directory>     Directory in which to find the images." << endl;
	cerr << "<-gaussian mean sigma>  Use Gaussian with mean and sigma for kernel-based smoothing\n";
	cerr << "<-epsilon value>        Epsilon value for ignoring image in kernel-based smoothing " << endl;
	cerr << "                          (default = " << EPSILON << ")" << endl;
	cerr << "<-scaling value>        Scaling value for final atlas (default = 1)\n";
	cerr << "<-norm>                 Normalise intensities before atlas construction " << endl;
	cerr << "<-maxnumber value>      Maxium number of images taken " << endl;
	cerr << "                          (default = off)\n";
	exit(1);
}

template<typename VoxelType> void normalise(irtkGenericImage<VoxelType> *input, irtkGenericImage<double> *output)
{
	int j;
	double mean, std;
	VoxelType *ptr1;
	double *ptr2;

	// Set up input histogram
	cout << "Setting up input histogram...";
	cout.flush();
	VoxelType min, max;
	input->GetMinMax(&min, &max);
	irtkHistogram histogram(max - min + 1);
	histogram.PutMin(min - 0.5);
	histogram.PutMax(max + 0.5);

	ptr1 = input->GetPointerToVoxels();
	for (j = 0; j < input->GetNumberOfVoxels(); j++) {
		histogram.AddSample(*ptr1);
		ptr1++;
	}

	mean = histogram.Mean();
	std = histogram.StandardDeviation();
	cout << "done" << endl;
	cout << "Stats: Mean = " << mean << " SD = " << std << endl;

	ptr1 = input->GetPointerToVoxels();
	ptr2 = output->GetPointerToVoxels();
	for (j = 0; j < input->GetNumberOfVoxels(); j++) {
		*ptr2 = (*ptr1 - mean) / std;
		ptr1++;
		ptr2++;
	}
}

void finalize_atlas(irtkGenericImage<double> *output, double scale, irtkGenericImage<double> *weight)
{
	int i;
	double *ptr2,*ptr3;

	cout << "Constructing average atlas...";
	cout.flush();
	ptr2 = output->GetPointerToVoxels();
	ptr3 = weight->GetPointerToVoxels();
	for (i = 0; i < output->GetNumberOfVoxels(); i++) {
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

template<typename VoxelType> void add_atlas(irtkGenericImage<VoxelType> *input, irtkGenericImage<double> *output, int norm, double inputweight, irtkGenericImage<double> *weight)
{
	int i;
	double *ptr1, *ptr2, *ptr3;

	cout << inputweight ;
	cout << " Adding input to atlas...";
	cout.flush();

	irtkGenericImage<double> tmp = *input;
	if (norm == true) normalise(input, &tmp);
	tmp *= inputweight;

	if (!(input->GetImageAttributes() == output->GetImageAttributes())
		||!(input->GetImageAttributes() == weight->GetImageAttributes())){
			cerr << "Mismatch of image attributes:" << endl;
			input->GetImageAttributes().Print();
			output->GetImageAttributes().Print();
			weight->GetImageAttributes().Print();
			exit(1);
	}

	ptr1 = tmp.GetPointerToVoxels();
	ptr2 = output->GetPointerToVoxels();
	ptr3 = weight->GetPointerToVoxels();

	for (i = 0; i < input->GetNumberOfVoxels(); i++) {
		if(*ptr1 >= 0){
			*ptr2 += *ptr1;
			*ptr3 += inputweight;
		}

		ptr1++;
		ptr2++;
		ptr3++;
	}

	cout << "done\n";
}

int main(int argc, char **argv)
{
	irtkImage                *input = NULL;
	irtkGenericImage<double> *output = NULL;
	irtkGenericImage<double> *weight = NULL;
	double scale, mean, sigma, epsilon;
	int padding, norm, i, no, ok, image_type, useGaussianWeights, maxiumimage = 10000;

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

	// Default: Use weights directly. Alternatively, convert using Gaussian formula.
	useGaussianWeights = false;

	// Default values for kernel smoothing
	mean = 0;
	sigma = 1;

	// Default: No intensity  normalisation
	norm = false;

	// Default: Epsilon
	epsilon = EPSILON;

	// Default padding
	padding = MIN_GREY;

	// Parse input file names and values
	char **input_name = new char *[MAX_IMAGES];
	double *input_value = new double[MAX_IMAGES];
	double *input_weight = new double[MAX_IMAGES];

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
		if ((ok == false) && (strcmp(argv[1], "-gaussian") == 0)) {
			argc--;
			argv++;
			useGaussianWeights = true;
			mean = atof(argv[1]);
			argc--;
			argv++;
			sigma = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-norm") == 0)) {
			argc--;
			argv++;
			norm = true;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-imagenames") == 0)) {
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
		if ((ok == false) && (strcmp(argv[1], "-padding") == 0)) {
			argc--;
			argv++;
			padding = atoi(argv[1]);
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
		if ((ok == false) && (strcmp(argv[1], "-epsilon") == 0)) {
			argc--;
			argv++;
			epsilon = atof(argv[1]);
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

					if (useGaussianWeights == true) {
						// Convert input value to weight based on Gaussian centred at above specified mean.
						input_weight[no] = 1.0 / sqrt(2.0) / sigma * exp(-pow((mean - input_value[no]) / sigma, 2.0) / 2);
					} else {
						// Use input value directly as a weight.
						input_weight[no] = input_value[no];
					}

					if (input_weight[no] > epsilon) {
						no++;
					}
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

		// Check what voxel type the first input image has: This voxel type will be used for all calculations
		if (i == 0) {
			irtkFileToImage *reader = irtkFileToImage::New(input_name[i]);
			image_type = reader->GetDataType();
		}

		switch (image_type) {
		case IRTK_VOXEL_UNSIGNED_CHAR: {
			cout << "Reading input image " << input_name[i] << "... ";
			cout.flush();
			input = new irtkGenericImage<unsigned char>(input_name[i]);
			cout << "done" << endl;

			if (i == 0) {
				output = new irtkGenericImage<double>(input->GetImageAttributes());
				weight = new irtkGenericImage<double>(input->GetImageAttributes());
			}

			add_atlas(dynamic_cast<irtkGenericImage<unsigned char> *>(input), output, norm, input_weight[i], weight);

			delete input;

			break;
									   }
		case IRTK_VOXEL_CHAR: {
			cout << "Reading input image " << input_name[i] << "... ";
			cout.flush();
			input = new irtkGenericImage<char>(input_name[i]);
			cout << "done" << endl;

			if (i == 0) {
				output = new irtkGenericImage<double>(input->GetImageAttributes());
				weight = new irtkGenericImage<double>(input->GetImageAttributes());
			}

			add_atlas(dynamic_cast<irtkGenericImage<char> *>(input), output, norm, input_weight[i], weight);

			delete input;

			break;
							  }
		case IRTK_VOXEL_UNSIGNED_SHORT: {
			cout << "Reading input image " << input_name[i] << "... ";
			cout.flush();
			input = new irtkGenericImage<unsigned short>(input_name[i]);
			cout << "done" << endl;

			if (i == 0) {
				output = new irtkGenericImage<double>(input->GetImageAttributes());
				weight = new irtkGenericImage<double>(input->GetImageAttributes());
			}

			add_atlas(dynamic_cast<irtkGenericImage<unsigned short> *>(input), output, norm, input_weight[i], weight);

			delete input;

			break;
										}
		case IRTK_VOXEL_SHORT: {
			cout << "Reading input image " << input_name[i] << "... ";
			cout.flush();
			input = new irtkGenericImage<short>(input_name[i]);
			cout << "done" << endl;

			if (i == 0) {
				output = new irtkGenericImage<double>(input->GetImageAttributes());
				weight = new irtkGenericImage<double>(input->GetImageAttributes());
			}

			add_atlas(dynamic_cast<irtkGenericImage<short> *>(input), output, norm, input_weight[i], weight);

			delete input;

			break;
							   }
		case IRTK_VOXEL_FLOAT: {
			cout << "Reading input image " << input_name[i] << "... ";
			cout.flush();
			input = new irtkGenericImage<float>(input_name[i]);
			cout << "done" << endl;

			if (i == 0) {
				output = new irtkGenericImage<double>(input->GetImageAttributes());
				weight = new irtkGenericImage<double>(input->GetImageAttributes());
			}

			add_atlas(dynamic_cast<irtkGenericImage<float> *>(input), output, norm, input_weight[i], weight);

			delete input;

			break;
							   }
		case IRTK_VOXEL_DOUBLE: {
			cout << "Reading input image " << input_name[i] << "... ";
			cout.flush();
			input = new irtkGenericImage<double>(input_name[i]);
			cout << "done" << endl;

			if (i == 0) {
				output = new irtkGenericImage<double>(input->GetImageAttributes());
				weight = new irtkGenericImage<double>(input->GetImageAttributes());
			}

			add_atlas(dynamic_cast<irtkGenericImage<double> *>(input), output, norm, input_weight[i], weight);

			delete input;

			break;
								}
		default:
			cerr << "Voxel type not supported" << endl;
			exit(1);
		}
	}

	// Finalize atlas
	finalize_atlas(output, scale, weight);

	// Write atlas
	cout << "Writing atlas to " << output_name << " ... ";
	cout.flush();
	output->Write(output_name);
	cout << "done\n";
}
