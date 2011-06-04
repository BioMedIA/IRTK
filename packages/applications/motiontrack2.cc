/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: motiontrack.cc 74 2009-08-13 11:38:17Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2009-08-13 12:38:17 +0100 (�? 13 八月 2009) $
  Version   : $Revision: 74 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkImage.h>

#include <irtkRegistration2.h>

#include <irtkGaussianBlurring.h>

#include <irtkGaussianBlurring4D.h>

char *dofout_name = NULL, *parin_name  = NULL, *parout_name = NULL, **filenames = NULL, *mask_name = NULL;;

void usage()
{
  cerr << "Usage: motiontrack [image sequence] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parin file>        Read parameter from file" << endl;
  cerr << "<-parout file>       Write parameter to file" << endl;
  cerr << "<-dofout file>       Write transformation to file" << endl;
  cerr << "<-ref file>          Reference time frame (default = first frame of image sequence)" << endl;
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
  exit(1);
}


int main(int argc, char **argv)
{
  int i, n, t, x, y, z, x1, y1, z1, t1, x2, y2, z2, t2, ok, debug;
  double spacing, sigma, xaxis[3], yaxis[3], zaxis[3];
  irtkGreyPixel padding;
  irtkImageFreeFormRegistrationMode mode;
  irtkMultiLevelFreeFormTransformation *mffd;

  // Check command line
  if (argc < 2) {
    usage();
  }

  // Get image names for sequence
  n = 0;
  filenames = argv;
  while ((argc > 1) && (argv[1][0] != '-' )) {
    argv++;
    argc--;
    n++;
  }
  filenames++;

  // Read image sequence
  cout << "Reading image sequence ... "; cout.flush();
  irtkGreyImage **image = new irtkGreyImage *[n];
  for (i = 0; i < n; i++) {
    cout << filenames[i] << endl;
    image[i] = new irtkGreyImage(filenames[i]);
  }

  // Check if there is at least one image
  if (n == 0) {
    usage();
  }

  // Fix ROI
  x1 = 0;
  y1 = 0;
  z1 = 0;
  t1 = 0;
  x2 = image[0]->GetX();
  y2 = image[0]->GetY();
  z2 = image[0]->GetZ();
  t2 = image[0]->GetT();

  // Default parameters
  padding   = MIN_GREY;
  spacing   = 0;
  sigma     = 0;
  mode      = RegisterXYZ;
  debug     = false;

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
  image[0]->GetOrientation(xaxis, yaxis, zaxis);

  // If there is an region of interest, use it
  if ((x1 != 0) || (x2 != image[0]->GetX()) ||
	  (y1 != 0) || (y2 != image[0]->GetY()) ||
	  (z1 != 0) || (z2 != image[0]->GetZ()) ||
	  (t1 != 0) || (t2 != image[0]->GetT())) {
		  for (i = 0; i < n; i++) {
			  *image[i] = image[i]->GetRegion(x1, y1, z1, t1, x2, y2, z2, t2);
		  }
  }

  // If sigma is larger than 0, blur images using 4D blurring
  if (sigma > 0) {
	  cout << "Blurring image sequences ... "; cout.flush();
	  for (i = 0; i < n; i++) {
		  irtkGaussianBlurring4D<irtkGreyPixel> gaussianBlurring4D(sigma);
		  gaussianBlurring4D.SetInput (image[i]);
		  gaussianBlurring4D.SetOutput(image[i]);
		  gaussianBlurring4D.Run();
	  }
	  cout << "done" << endl;
  }

  // Use identity transformation to start
  mffd = new irtkMultiLevelFreeFormTransformation;

  for( i = 0; i < n; i ++){
	  for (t = 1; t < image[i]->GetT(); t++) {

		  // Create registration filter
		  irtkImageFreeFormRegistration2 *registration = NULL;
		  if (image[i]->GetZ() == 1) {
			  cerr<<"no 2D for motiontrack2 now"<<endl;
			  exit(1);
		  } else {
			  registration = new irtkImageFreeFormRegistration2;
		  }

		  // Combine images
		  irtkImageAttributes attr = image[i]->GetImageAttributes();
		  attr._t = 1;
		  irtkGreyImage *target = new irtkGreyImage(attr);
		  irtkGreyImage *source = new irtkGreyImage(attr);
		  for (z = 0; z < target->GetZ(); z++) {
			  for (y = 0; y < target->GetY(); y++) {
				  for (x = 0; x < target->GetX(); x++) {
					  target->Put(x, y, z, 0, image[i]->Get(x, y, z, 0));
					  source->Put(x, y, z, 0, image[i]->Get(x, y, z, t));
				  }
			  }
		  }

          // Mask target
          if (mask_name != NULL) {
              irtkGreyImage mask(mask_name);
              irtkInterpolateImageFunction *interpolator = new irtkShapeBasedInterpolateImageFunction;
              double xsize, ysize, zsize, size;
              mask.GetPixelSize(&xsize, &ysize, &zsize);
              size = xsize;
              size = (size < ysize) ? size : ysize;
              size = (size < zsize) ? size : zsize;
              irtkResampling<irtkGreyPixel> resampling(size, size, size);
              resampling.SetInput (&mask);
              resampling.SetOutput(&mask);
              resampling.SetInterpolator(interpolator);
              resampling.Run();
              delete interpolator;
              for (z = 0; z < target->GetZ(); z++) {
                  for (y = 0; y < target->GetY(); y++) {
                      for (x = 0; x < target->GetX(); x++) {
                          double wx,wy,wz;
                          wx = x; wy = y; wz = z;
                          target->ImageToWorld(wx,wy,wz);
                          mask.WorldToImage(wx,wy,wz);
                          wx = round(wx); wy = round(wy); wz = round(wz);
                          if(wx >= 0 && wx < mask.GetX() && wy >= 0 && wy < mask.GetY()
                              && wz >= 0 && wz < mask.GetZ() && mask.GetAsDouble(wx,wy,wz) <= 0){
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

		  // Read parameter if there any, otherwise make an intelligent guess
		  if (parin_name != NULL) {
			  registration->irtkImageRegistration2::Read(parin_name);
		  } else {
			  registration->GuessParameter();
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

		  // Write parameters if necessary
		  if (parout_name != NULL) {
			  registration->irtkImageRegistration2::Write(parout_name);
		  }

		  // Run registration filter
		  registration->Run();

		  if(registration->GetMFFDMode()){
			  mffd->CombineLocalTransformation();
		  }

		  // Write the final transformation estimate
		  if (dofout_name != NULL) {
			  char buffer[255];
			  sprintf(buffer, "%s\\%d_sequence_%.2d.dof.gz", dofout_name,i, t);
			  mffd->irtkTransformation::Write(buffer);
		  }

		  delete registration;
		  delete target;
		  delete source;
	  }
  }
}
