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

#include <irtkImageFunction.h>

#include <irtkHistogram.h>

#include <irtkTransformation.h>

// Default filenames
char *source_name = NULL, *target_name = NULL, *mask_name = NULL;
char *trans_name  = NULL, *output_name = NULL;

#define DEFAULT_BINS 255

// Default number of bins
int nbins_x = 0, nbins_y = 0;

void usage()
{
  cerr << "Usage: evaluation [target] [source] <options>\n" << endl;
  cerr << "where <options> is one or more of the following: \n" << endl;
  cerr << "<-dofin file>      Input transformation" << endl;
  cerr << "<-nbins_x no>      Number of bins for target (Default smaller of dynamic range or 255)" << endl;
  cerr << "<-nbins_y no>      Number of bins for source (Default as above)" << endl;
  cerr << "<-Rx1 pixel>       Region of interest" << endl;
  cerr << "<-Ry1 pixel>       Region of interest" << endl;
  cerr << "<-Rz1 pixel>       Region of interest" << endl;
  cerr << "<-Rx2 pixel>       Region of interest" << endl;
  cerr << "<-Ry2 pixel>       Region of interest" << endl;
  cerr << "<-Rz2 pixel>       Region of interest" << endl;
  cerr << "<-Tp  value>       Padding value in target" << endl;
  cerr << "<-mask file>       Binary mask to define ROI" << endl;
  cerr << "<-linear>          Linear interpolation" << endl;
  cerr << "<-bspline>         B-spline interpolation" << endl;
  cerr << "<-cspline>         Cubic spline interpolation" << endl;
  cerr << "<-sinc>            Sinc interpolation" << endl;
  cerr << "<-output>          Output similarity values to file" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation *transformation = NULL;
  irtkInterpolateImageFunction *interpolator = NULL;
  irtkGreyPixel target_min, source_min, target_max, source_max;
  int ok, i, x, y, z, i1, j1, k1, i2, j2, k2;
  double x1, y1, z1, x2, y2, z2, Tp, widthx, widthy, val;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  // Read target and source image
  irtkGreyImage target(target_name);
  irtkGreyImage source(source_name);

  // Default padding
  Tp = -1.0 * FLT_MAX;

  // Fix ROI
  i1 = 0;
  j1 = 0;
  k1 = 0;
  i2 = target.GetX();
  j2 = target.GetY();
  k2 = target.GetZ();

  // Fix no. of bins;
  nbins_x = 0;
  nbins_y = 0;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      trans_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nbins_x") == 0)) {
      argc--;
      argv++;
      nbins_x = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nbins_y") == 0)) {
      argc--;
      argv++;
      nbins_y = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      Tp = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      i1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      i2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      j1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      j2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      k1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      k2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-linear") == 0)) {
      argc--;
      argv++;
      if (source.GetZ() == 1){
        interpolator = new irtkLinearInterpolateImageFunction2D;
      } else {
        interpolator = new irtkLinearInterpolateImageFunction;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)) {
      argc--;
      argv++;
      if (source.GetZ() == 1){
      	interpolator = new irtkBSplineInterpolateImageFunction2D;
      } else {
      	interpolator = new irtkBSplineInterpolateImageFunction;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-cspline") == 0)) {
      argc--;
      argv++;
      if (source.GetZ() == 1){
      	interpolator = new irtkCSplineInterpolateImageFunction2D;
      } else {
      	interpolator = new irtkCSplineInterpolateImageFunction;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sinc") == 0)) {
      argc--;
      argv++;
      if (source.GetZ() == 1){
      	interpolator = new irtkSincInterpolateImageFunction2D;
      } else {
      	interpolator = new irtkSincInterpolateImageFunction;
      }
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
    if ((ok == false) && (strcmp(argv[1], "-output") == 0)) {
        argc--;
        argv++;
        output_name = argv[1];
        argc--;
        argv++;
        ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkGreyImage mask;

  // Set up a mask,
  if (mask_name == NULL) {
    mask.Read(target_name);

    irtkGreyPixel *ptr2mask = mask.GetPointerToVoxels();
    irtkGreyPixel *ptr2tgt  = target.GetPointerToVoxels();

    for (i = 0; i < target.GetNumberOfVoxels(); i++) {
      if (*ptr2tgt > Tp)
        *ptr2mask = 1;
      else
        *ptr2mask = 0;

      ++ptr2tgt;
      ++ptr2mask;
    }
  } else {
    mask.Read(mask_name);
    if ((mask.GetX() != target.GetX()) ||
        (mask.GetY() != target.GetY()) ||
        (mask.GetZ() != target.GetZ())) {
      cerr << "evaluation2: Target and mask dimensions mismatch. Exiting." << endl;
      exit(1);
    }
  }

  // If there is an region of interest, use it
  if ((i1 != 0) || (i2 != target.GetX()) ||
      (j1 != 0) || (j2 != target.GetY()) ||
      (k1 != 0) || (k2 != target.GetZ())) {
    target = target.GetRegion(i1, j1, k1, i2, j2, k2);
    source = source.GetRegion(i1, j1, k1, i2, j2, k2);
    mask   = mask.GetRegion(i1, j1, k1, i2, j2, k2);
  }

  // Set min and max of histogram
  target.GetMinMax(&target_min, &target_max);
  source.GetMinMax(&source_min, &source_max);

  // Calculate number of bins to use
  if (nbins_x == 0) {
    nbins_x = (int) round(target_max - target_min) + 1;
    if (nbins_x > DEFAULT_BINS)
      nbins_x = DEFAULT_BINS;
  }

  if (nbins_y == 0) {
    nbins_y = (int) round(source_max - source_min) + 1;
    if (nbins_y > DEFAULT_BINS)
      nbins_y = DEFAULT_BINS;
  }

  // Create default interpolator if necessary
  if (interpolator == NULL) {
  	if (source.GetZ() == 1){
  		interpolator = new irtkNearestNeighborInterpolateImageFunction2D;
  	} else {
  		interpolator = new irtkNearestNeighborInterpolateImageFunction;
  	}
  }
  interpolator->SetInput(&source);
  interpolator->Initialize();

  // Calculate the source image domain in which we can interpolate
  interpolator->Inside(x1, y1, z1, x2, y2, z2);

  // Create histogram
  irtkHistogram_2D<int> histogram(nbins_x, nbins_y);
  widthx = (target_max - target_min) / (nbins_x - 1.0);
  widthy = (source_max - source_min) / (nbins_y - 1.0);

  histogram.PutMin(target_min - 0.5*widthx, source_min - 0.5*widthy);
  histogram.PutMax(target_max + 0.5*widthx, source_max + 0.5*widthy);

  if (trans_name == 0) {
    transformation = new irtkRigidTransformation;
  } else {
    transformation = irtkTransformation::New(trans_name);
  }

  target_min = FLT_MAX;
  source_min = FLT_MAX;
  target_max = -1.0 * FLT_MAX;
  source_max = -1.0 * FLT_MAX;

  // Fill histogram
  for (z = 0; z < target.GetZ(); z++) {
    for (y = 0; y < target.GetY(); y++) {
      for (x = 0; x < target.GetX(); x++) {

        if (mask(x, y, z) > 0) {
          val = target(x, y, z);

          if (val > target_max)
            target_max = val;
          if (val < target_min)
            target_min = val;

          irtkPoint p(x, y, z);
          // Transform point into world coordinates
          target.ImageToWorld(p);
          // Transform point
          transformation->Transform(p);
          // Transform point into image coordinates
          source.WorldToImage(p);

        	// A bad thing might happen for the 2D case.
        	if ((source.GetZ() == 1) &&
		    (p._z > 0.5 || p._z < -0.5)){
		  cerr << "Transformed point outside plane of 2D source image." << endl;
		  exit(1);
        	}

        	// 2D and in plane but out of FoV.
        	if ((source.GetZ() == 1) &&
		    (p._x <= x1 || p._x >= x2 ||
		     p._y <= y1 || p._y >= y2))
		  continue;

        	// 3D and out of FoV.
        	if ((source.GetZ() > 1) &&
		    (p._x <= x1 || p._x >= x2 ||
		     p._y <= y1 || p._y >= y2 ||
		     p._z <= z1 || p._z >= z2))
		  continue;

		// Should be able to interpolate if we've got this far.

        	val = interpolator->EvaluateInside(p._x, p._y, p._z);

			if(val > Tp){
        		histogram.AddSample(target(x, y, z), val);
			}

        	if (val >  source_max)
        		source_max = val;
        	if (val < source_min)
        		source_min = val;

        }
      }
    }
  }

  cout << "SSD: "    << histogram.SumsOfSquaredDifferences() / (double)histogram.NumberOfSamples() << endl;
  cout << "CC: "     << histogram.CrossCorrelation() << endl;
  cout << "NMI: "    << histogram.NormalizedMutualInformation() << endl;
  if (nbins_x == nbins_y) {
      cout << "Label consistency: " << histogram.LabelConsistency() << endl;
      cout << "Kappa statistic: " << histogram.Kappa() << endl;
    }

  if(output_name){
      cerr << "Writing Results of SSD CC and NMI to " << output_name << endl;
      ofstream fout(output_name,ios::app);
      fout << histogram.SumsOfSquaredDifferences() / (double)histogram.NumberOfSamples() <<endl;
      fout << histogram.CrossCorrelation() <<endl;
      fout << histogram.NormalizedMutualInformation() <<endl;
	  if (nbins_x == nbins_y) {
		  fout << histogram.LabelConsistency() << endl;
	  }
      fout.close();
  }
}
