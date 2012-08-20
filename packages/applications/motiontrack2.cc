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

#include <irtkRegistration2.h>

#include <irtkGaussianBlurring.h>

#include <irtkGaussianBlurring4D.h>

char *dofout_name = NULL, *parin_name = NULL, *parout_name = NULL, *image_name = NULL, *mask_name =
		NULL;
;

void usage()
{
	cerr << "Usage: motiontrack [image sequence] <options> \n" << endl;
	cerr << "where <options> is one or more of the following:\n" << endl;
	cerr << "<-parin file>        Read parameter from file" << endl;
	cerr << "<-parout file>       Write parameter to file" << endl;
	cerr << "<-dofout file>       Write transformation to file" << endl;
	cerr << "<-ref file>          Reference time frame (default = first frame of image sequence)"
			<< endl;
	cerr << "<-Rx1 value>         Region of interest in images" << endl;
	cerr << "<-Ry1 value>         Region of interest in images" << endl;
	cerr << "<-Rz1 value>         Region of interest in images" << endl;
	cerr << "<-Rt1 value>         Region of interest in images" << endl;
	cerr << "<-Rx2 value>         Region of interest in images" << endl;
	cerr << "<-Ry2 value>         Region of interest in images" << endl;
	cerr << "<-Rz2 value>         Region of interest in images" << endl;
	cerr << "<-Rt2 value>         Region of interest in images" << endl;
	cerr << "<-Tp  value>         Padding value" << endl;
	cerr << "<-mask file>         Use a mask to define the ROI. The mask" << endl;
	cerr << "<-adaptive weight>   Adapt regulation weight using the weight" << endl;
	cerr << "<-mode 0/1>          Registration mode [t0-ti]/[ti-ti+1]" << endl;
	cerr << "                     (w = weight^(t (t<n/2) n-t (t>-n/2))*w)" << endl;
	exit(1);
}

int main(int argc, char **argv)
{
	int t, x, y, z, x1, y1, z1, t1, x2, y2, z2, t2, framemode, ok, debug;
	double spacing, sigma, adaptive, xaxis[3], yaxis[3], zaxis[3];
	irtkGreyPixel padding;
	irtkImageFreeFormRegistrationMode mode;
	irtkMultiLevelFreeFormTransformation *mffd;

	// Check command line
	if (argc < 2) {
		usage();
	}

	// Get image names for sequence
	image_name = argv[1];
	argv++;
	argc--;

	// Read image sequence
	cout << "Reading image sequence ... ";
	cout.flush();
	irtkGreyImage *image = new irtkGreyImage(image_name);

	// Fix ROI
	x1 = 0;
	y1 = 0;
	z1 = 0;
	t1 = 0;
	x2 = image->GetX();
	y2 = image->GetY();
	z2 = image->GetZ();
	t2 = image->GetT();

	// Default parameters
	padding = MIN_GREY;
	spacing = 0;
	sigma = 0;
	mode = RegisterXYZ;
	debug = false;
	adaptive = 0;
	framemode = 0;

	// Parse remaining parameters
	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
			argc--;
			argv++;
			x1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
			argc--;
			argv++;
			x2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
			argc--;
			argv++;
			y1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
			argc--;
			argv++;
			y2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
			argc--;
			argv++;
			z1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
			argc--;
			argv++;
			z2 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rt1") == 0)) {
			argc--;
			argv++;
			t1 = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rt2") == 0)) {
			argc--;
			argv++;
			t2 = atoi(argv[1]);
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

		if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
			argc--;
			argv++;
			padding = atoi(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-mask") == 0)) {
			argc--;
			argv++;
			mask_name = argv[1];
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-ds") == 0)) {
			argc--;
			argv++;
			spacing = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-parin") == 0)) {
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
		if ((ok == false) && (strcmp(argv[1], "-blur") == 0)) {
			argc--;
			argv++;
			ok = true;
			sigma = atof(argv[1]);
			argc--;
			argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-adaptive") == 0)) {
			argc--;
			argv++;
			ok = true;
			adaptive = atof(argv[1]);
			argc--;
			argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-mode") == 0)) {
			argc--;
			argv++;
			ok = true;
			framemode = atoi(argv[1]);
			argc--;
			argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-debug") == 0)) {
			argc--;
			argv++;
			ok = true;
			debug = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-xy_only") == 0)) {
			argc--;
			argv++;
			mode = RegisterXY;
			ok = true;
		}

		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	// Image orientation
	image->GetOrientation(xaxis, yaxis, zaxis);

	// If there is an region of interest, use it
	if ((x1 != 0) || (x2 != image->GetX()) ||
			(y1 != 0) || (y2 != image->GetY()) ||
			(z1 != 0) || (z2 != image->GetZ()) ||
			(t1 != 0) || (t2 != image->GetT())) {
		*image = image->GetRegion(x1, y1, z1, t1, x2, y2, z2, t2);
	}

	// If sigma is larger than 0, blur images using 4D blurring
	if (sigma > 0) {
		cout << "Blurring image sequences ... ";
		cout.flush();
		irtkGaussianBlurring4D<irtkGreyPixel> gaussianBlurring4D(sigma);
		gaussianBlurring4D.SetInput(image);
		gaussianBlurring4D.SetOutput(image);
		gaussianBlurring4D.Run();
		cout << "done" << endl;
	}

	// Use identity transformation to start
	mffd = new irtkMultiLevelFreeFormTransformation;

	for (t = 1; t < image->GetT(); t++) {

		// Create registration filter
		irtkImageFreeFormRegistration2 *registration = new irtkImageFreeFormRegistration2;

		// if frame to frame clear previous mffd
		if (framemode != 0) {
			delete mffd;
			mffd = new irtkMultiLevelFreeFormTransformation;
		}

		// Combine images
		irtkImageAttributes attr = image->GetImageAttributes();
		attr._t = 1;
		irtkGreyImage *target = new irtkGreyImage(attr);
		irtkGreyImage *source = new irtkGreyImage(attr);
		for (z = 0; z < target->GetZ(); z++) {
			for (y = 0; y < target->GetY(); y++) {
				for (x = 0; x < target->GetX(); x++) {
					if (framemode == 0) {
						target->Put(x, y, z, 0, image->Get(x, y, z, 0));
					}
					else {
						target->Put(x, y, z, 0, image->Get(x, y, z, t - 1));
					}
					source->Put(x, y, z, 0, image->Get(x, y, z, t));
				}
			}
		}

		// Mask target
		irtkGreyImage mask;
		if (mask_name != NULL) {
			mask.Read(mask_name);
			irtkInterpolateImageFunction *interpolator = new irtkShapeBasedInterpolateImageFunction;
			double xsize, ysize, zsize, size;
			mask.GetPixelSize(&xsize, &ysize, &zsize);
			size = xsize;
			size = (size < ysize) ? size : ysize;
			size = (size < zsize) ? size : zsize;
			irtkResampling<irtkGreyPixel> resampling(size, size, size);
			resampling.SetInput(&mask);
			resampling.SetOutput(&mask);
			resampling.SetInterpolator(interpolator);
			resampling.Run();
			delete interpolator;
			for (z = 0; z < target->GetZ(); z++) {
				for (y = 0; y < target->GetY(); y++) {
					for (x = 0; x < target->GetX(); x++) {
						double wx, wy, wz;
						wx = x;
						wy = y;
						wz = z;
						target->ImageToWorld(wx, wy, wz);
						mask.WorldToImage(wx, wy, wz);
						wx = round(wx);
						wy = round(wy);
						wz = round(wz);
						if (wx >= 0 && wx < mask.GetX() && wy >= 0 && wy < mask.GetY()
								&& wz >= 0 && wz < mask.GetZ() && mask.GetAsDouble(wx, wy, wz) <= 0) {
							target->Put(x, y, z, 0, padding);
						}
					}
				}
			}
		}

		// Set input and output for the registration filter
		registration->SetInput(target, source);
		registration->SetOutput(mffd);
		registration->SetMode(mode);
		registration->SetDebugFlag(debug);

		// Make an initial Guess for the parameters.
		registration->GuessParameter();
		// Overrride with any the user has set.
		if (parin_name != NULL) {
			registration->irtkImageRegistration2::Read(parin_name);
		}

		// Override parameter settings if necessary
		if (padding != MIN_GREY) {
			registration->SetTargetPadding(padding);
		}
		if (spacing > 0) {
			registration->SetDX(spacing);
			registration->SetDY(spacing);
			registration->SetDZ(spacing);
		}
		if (adaptive > 0) {
			int times = 0;
			t < round(image->GetT() / 3) ? (times = t) : (times = round(image->GetT() / 3));
			if (t > round(image->GetT() * 2 / 3)) times = (image->GetT() - 1 - t);
			double weight = registration->GetLambda1();
			weight = pow(adaptive, times) * weight;
			cout << "current lambda1 of frame " << t << " is " << weight << endl;
			registration->SetLambda1(weight);
		}

		// Write parameters if necessary
		if (parout_name != NULL) {
			registration->irtkImageRegistration2::Write(parout_name);
		}

		// Run registration filter
		registration->Run();

		// Write the final transformation estimate
		if (dofout_name != NULL) {
			char buffer[255];
			sprintf(buffer, "%s%.2d.dof.gz", dofout_name, t);
			if (framemode == 0) {
				if (registration->GetMFFDMode()) {
					mffd->CombineLocalTransformation();
				}
				mffd->irtkTransformation::Write(buffer);
			}
			else {
				mffd->irtkTransformation::Write(buffer);
			}
		}
		else {
			if (registration->GetMFFDMode()) {
				mffd->CombineLocalTransformation();
			}
		}

		delete registration;
		delete target;
		delete source;
	}
}
