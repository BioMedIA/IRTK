/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkRegistration2.h>

void usage()
{
	cerr << "Usage: computeTFFD [target] [N] [source] [t_real] <options> \n" << endl;
	cerr << "takes in a reference image [source] and [N] [target] images with \n" << endl;
	cerr << "a text file [t_real] containing the time points of the target image\n" << endl;
	cerr << "whith <options> is one or more of the following:\n" << endl;
	cerr << "<-parin file>        Read parameter from file" << endl;
	cerr << "<-parout file>       Write parameter to file" << endl;
	cerr << "<-dofin  file>       Read transformation from file" << endl;
	cerr << "<-dofout file>       Write transformation to file" << endl;
	cerr << "<-Rx1 value>         Region of interest in both images" << endl;
	cerr << "<-Ry1 value>         Region of interest in both images" << endl;
	cerr << "<-Rz1 value>         Region of interest in both images" << endl;
	cerr << "<-Rx2 value>         Region of interest in both images" << endl;
	cerr << "<-Ry2 value>         Region of interest in both images" << endl;
	cerr << "<-Rz2 value>         Region of interest in both images" << endl;
	cerr << "<-Tx1 value>         Region of interest in target image" << endl;
	cerr << "<-Ty1 value>         Region of interest in target image" << endl;
	cerr << "<-Tz1 value>         Region of interest in target image" << endl;
	cerr << "<-Tx2 value>         Region of interest in target image" << endl;
	cerr << "<-Ty2 value>         Region of interest in target image" << endl;
	cerr << "<-Tz2 value>         Region of interest in target image" << endl;
	cerr << "<-Sx1 value>         Region of interest in source image" << endl;
	cerr << "<-Sy1 value>         Region of interest in source image" << endl;
	cerr << "<-Sz1 value>         Region of interest in source image" << endl;
	cerr << "<-Sx2 value>         Region of interest in source image" << endl;
	cerr << "<-Sy2 value>         Region of interest in source image" << endl;
	cerr << "<-Sz2 value>         Region of interest in source image" << endl;
	cerr << "<-Tp  value>         Padding value in target image" << endl;
	cerr << "<-ds  value>         Initial control point spacing" << endl;
	cerr << "<-dsT value>         Initial temporal control point spacing" << endl;
	cerr << "<-Sp  value>         Smoothness preservation value" << endl;
	cerr << "<-Gn  value>         Gradient normalize value" << endl;
	cerr << "<-debug>             Enable debugging information" << endl;
	cerr << "<-mask file>         Use a mask to define the ROI. The mask" << endl;
	cerr << "                     must have the same dimensions as the target." << endl;
	cerr << "                     Voxels in the mask with zero or less are " << endl;
	cerr << "                     padded in the target." << endl;
	cerr << "<-mask_dilation n>   Dilate mask n times before using it" << endl;

	exit(1);
}

int main(int argc, char **argv)
{
	double spacing, spacingT;
	int ok, padding, mask_dilation = 0;
	int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
	int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;
	double sp,gn;
	int N_source;
	double *t_real = NULL;

	// Default filenames
	char *target_name = NULL, **source_name = NULL, *timeFile_name = NULL;
	char *dofin_name  = NULL, *dofout_name = NULL;
	char *parin_name  = NULL, *parout_name = NULL;
	char *mask_name = NULL;

	// Check command line
	if (argc < 5) {
		usage();
	}

	// Parse source and target images
	target_name = argv[1];
	argc--;
	argv++;
	N_source = atoi(argv[1]);
	argc--;
	argv++;
	source_name = new char *[N_source];
	for (int n = 0; n < N_source; n++) {
		source_name[n] = argv[1];
		argc--;
		argv++;
	}
	timeFile_name = argv[1];
	argc--;
	argv++;

	// Read source image
	cout << "Reading target ... "; cout.flush();
	irtkGreyImage target(target_name);
	cout << "done" << endl;
	// Read target image
	cout << "Reading source ... "; cout.flush();
	irtkGreyImage **source = new irtkGreyImage *[N_source];
	for (int n = 0; n < N_source; n++) {
		source[n] = new irtkGreyImage(source_name[n]);
	}
	cout << "done" << endl;
	// Read text file containing time points for targets
	cout << "Reading time text file ... "; cout.flush();
	t_real = new double [N_source];
	ifstream textFileReader;
	string STRING;
	textFileReader.open(timeFile_name);
	if (textFileReader.is_open()) {
		for (int n = 0; n < N_source; n++) {
			if (!textFileReader.eof()) {
				getline(textFileReader,STRING);
				char * line = new char[STRING.size()];
				for (int j=0; j<STRING.size(); j++) line[j] = STRING[j];
				t_real[n] = atof(line);
			}
		}
	} else {
		cerr<<"coudn't open time file "<<timeFile_name<<endl;
		exit(1);
	}
	cout << "done" << endl;

	// Fix ROI
	target_x1 = 0;
	target_y1 = 0;
	target_z1 = 0;
	target_x2 = target.GetX();
	target_y2 = target.GetY();
	target_z2 = target.GetZ();
	source_x1 = 0;
	source_y1 = 0;
	source_z1 = 0;
	source_x2 = source[0]->GetX();
	source_y2 = source[0]->GetY();
	source_z2 = source[0]->GetZ();
	sp = 0;
	gn = 0;

	// Create registration filter
	irtkImageTSFFDRegistration *registration = NULL;
	registration = new irtkImageTSFFDRegistration;

	// Create initial multi-level free-form deformation
	irtkMultiLevelFreeFormTransformation *mffd = NULL;

	// Default parameters
	padding   = MIN_GREY;
	spacing   = 0;
	spacingT   = 0;

	// Parse remaining parameters
	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
			argc--;
			argv++;
			dofin_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-dofout") == 0)) {
			argc--;
			argv++;
			dofout_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
			argc--;
			argv++;
			target_x1 = source_x1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
			argc--;
			argv++;
			target_x2 = source_x2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
			argc--;
			argv++;
			target_y1 = source_y1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
			argc--;
			argv++;
			target_y2 = source_y2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
			argc--;
			argv++;
			target_z1 = source_z1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
			argc--;
			argv++;
			target_z2 = source_z2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tx1") == 0)) {
			argc--;
			argv++;
			target_x1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tx2") == 0)) {
			argc--;
			argv++;
			target_x2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ty1") == 0)) {
			argc--;
			argv++;
			target_y1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ty2") == 0)) {
			argc--;
			argv++;
			target_y2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tz1") == 0)) {
			argc--;
			argv++;
			target_z1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tz2") == 0)) {
			argc--;
			argv++;
			target_z2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sx1") == 0)) {
			argc--;
			argv++;
			source_x1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sx2") == 0)) {
			argc--;
			argv++;
			source_x2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sy1") == 0)) {
			argc--;
			argv++;
			source_y1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sy2") == 0)) {
			argc--;
			argv++;
			source_y2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sz1") == 0)) {
			argc--;
			argv++;
			source_z1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Sz2") == 0)) {
			argc--;
			argv++;
			source_z2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
			argc--;
			argv++;
			padding = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-debug") == 0)) {
			argc--;
			argv++;
			ok = true;
			registration->SetDebugFlag(true);
		}
		if ((ok == false) && (strcmp(argv[1], "-Sp") == 0)) {
			argc--;
			argv++;
			sp = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Gn") == 0)) {
			argc--;
			argv++;
			gn = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-x_only") == 0)) {
			argc--;
			argv++;
			registration->SetMode(RegisterX);
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-y_only") == 0)) {
			argc--;
			argv++;
			registration->SetMode(RegisterY);
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-xy_only") == 0)) {
			argc--;
			argv++;
			registration->SetMode(RegisterXY);
			ok = true;
		}
		if ((ok == false) && ((strcmp(argv[1], "-parameter") == 0) || (strcmp(argv[1], "-parin") == 0))) {
			argc--;
			argv++;
			ok = true;
			parin_name = argv[1];
			argc--;
			argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-parout") == 0)) {
			argc--;
			argv++;
			ok = true;
			parout_name = argv[1];
			argc--;
			argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-mask") == 0)) {
			argc--;
			argv++;
			mask_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-mask_dilation") == 0)) {
			argc--;
			argv++;
			mask_dilation = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	// Is there a mask to use?
	if (mask_name != NULL) {
		double xw, yw, zw, xwm, ywm, zwm;
		irtkGreyImage mask(mask_name);

		if (target.GetX() == 1){
			if (mask.GetY() != target.GetY() ||
				mask.GetZ() != target.GetZ()) {
					cerr << "Mask given does not match target dimensions." << endl;
					exit(1);
			}
		}
		else if (target.GetY() == 1){
			if (mask.GetX() != target.GetX() ||
				mask.GetZ() != target.GetZ()) {
					cerr << "Mask given does not match target dimensions." << endl;
					exit(1);
			}
		}
		else if (target.GetZ() == 1){
			if (mask.GetX() != target.GetX() ||
				mask.GetY() != target.GetY()) {
					cerr << "Mask given does not match target dimensions." << endl;
					exit(1);
			}
		}

		if (mask_dilation > 0) {
			irtkDilation<irtkGreyPixel> dilation;
			dilation.SetConnectivity(CONNECTIVITY_26);
			dilation.SetInput(&mask);
			dilation.SetOutput(&mask);
			cout << "Dilating mask ... ";
			cout.flush();
			for (int i = 0; i < mask_dilation; i++) dilation.Run();
			cout << "done" << endl;

		}

		for (int x = 0; x < target.GetX(); x++) {
			for (int y = 0; y < target.GetY(); y++) {
				for (int z = 0; z < target.GetZ(); z++) {
					xw = x;
					yw = y;
					zw = z;
					target.ImageToWorld(xw, yw, zw);
					xwm = xw;
					ywm = yw;
					zwm = zw;
					mask.WorldToImage(xwm, ywm, zwm);
					if (mask.Get(xwm, ywm, zwm) <= 0) {
						target.Put(x, y, z, padding);
					}
				}
			}
		}
	}

	// If there is an region of interest for the target image, use it
	if ((target_x1 != 0) || (target_x2 != target.GetX()) ||
		(target_y1 != 0) || (target_y2 != target.GetY()) ||
		(target_z1 != 0) || (target_z2 != target.GetZ())) {
			target = target.GetRegion(target_x1, target_y1, target_z1,
				target_x2, target_y2, target_z2);
	}

	for (int n = 0; n < N_source; n++) {
		// If there is an region of interest for the source image, use it
		if ((source_x1 != 0) || (source_x2 != source[n]->GetX()) ||
			(source_y1 != 0) || (source_y2 != source[n]->GetY()) ||
			(source_z1 != 0) || (source_z2 != source[n]->GetZ())) {
				*source[n] = source[n]->GetRegion(source_x1, source_y1, source_z1,
					source_x2, source_y2, source_z2);
		}
	}

	if (dofin_name != NULL) {
		irtkTransformation *transform = irtkTransformation::New(dofin_name);
		if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
			mffd = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
		} else {
			if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
				mffd = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
			} else {
				if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
					mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
				} else {
					cerr << "Input transformation is not of type rigid or affine " << endl;
					cerr << "or multi-level free form deformation" << endl;
					exit(1);
				}
			}
		}
		delete transform;
	} else {
		// Otherwise use identity transformation to start
		mffd = new irtkMultiLevelFreeFormTransformation;
	}

	// Set input and output for the registration filter
	registration->SetInput(&target, source, N_source);
	registration->SetOutput(mffd);
	registration->SetTime(t_real);

	// Make an initial Guess for the parameters.
	registration->GuessParameter();
	// Overrride with any the user has set.
	if (parin_name != NULL) {
		registration->irtkTemporalImageRegistration::Read(parin_name);
	}

	if(sp != 0) registration->SetLambda1(sp);
	if(gn != 0) registration->SetLambda3(gn);

	// Override parameter settings if necessary
	if (padding != MIN_GREY) {
		registration->SetTargetPadding(padding);
	}

	// Write parameters if necessary
	if (parout_name != NULL) {
		registration->irtkTemporalImageRegistration::Write(parout_name);
	}

	// Run registration filter
	registration->Run();

	// Write the final transformation estimate
	if (dofout_name != NULL) {
		mffd->irtkTransformation::Write(dofout_name);
	}

	delete registration;
	delete mffd;
	delete [] source;
}