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
char *from_name = NULL, *target_name = NULL, *output_name = NULL;
char *trans_name  = NULL;

#define DEFAULT_BINS 64

// Default number of bins
int nbins_x = 0, nbins_y = 0;

void usage()
{
  cerr << "Usage: temporalextract [from] [target] [output] <options>\n" << endl;
  cerr << "where <options> is one or more of the following: \n" << endl;
  cerr << "<-dofin file>      Input transformation" << endl;
  cerr << "<-linear>          Linear interpolation" << endl;
  cerr << "<-bspline>         B-spline interpolation" << endl;
  cerr << "<-cspline>         Cubic spline interpolation" << endl;
  cerr << "<-sinc>            Sinc interpolation" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkTransformation *transformation = NULL;
  irtkInterpolateImageFunction *interpolator = NULL;
  irtkRealPixel target_min, target_max, source_min, source_max;
  int ok, x, y, z, t;
  double x1, y1, z1, x2, y2, z2, widthx, widthy, val;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse target and target images
  from_name = argv[1];
  argc--;
  argv++;
  target_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read target and target image
  irtkRealImage target(target_name);
  irtkRealImage from(from_name);

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
    if ((ok == false) && (strcmp(argv[1], "-linear") == 0)) {
      argc--;
      argv++;
      if (target.GetZ() == 1){
        interpolator = new irtkLinearInterpolateImageFunction2D;
      } else {
        interpolator = new irtkLinearInterpolateImageFunction;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)) {
      argc--;
      argv++;
      if (target.GetZ() == 1){
      	interpolator = new irtkBSplineInterpolateImageFunction2D;
      } else {
      	interpolator = new irtkBSplineInterpolateImageFunction;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-cspline") == 0)) {
      argc--;
      argv++;
      if (target.GetZ() == 1){
      	interpolator = new irtkCSplineInterpolateImageFunction2D;
      } else {
      	interpolator = new irtkCSplineInterpolateImageFunction;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sinc") == 0)) {
      argc--;
      argv++;
      if (target.GetZ() == 1){
      	interpolator = new irtkSincInterpolateImageFunction2D;
      } else {
      	interpolator = new irtkSincInterpolateImageFunction;
      }
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }


  // Set min and max of histogram
  from.GetMinMax(&target_min, &target_max);
  target.GetMinMax(&source_min, &source_max);

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
  	if (target.GetZ() == 1){
  		interpolator = new irtkNearestNeighborInterpolateImageFunction2D;
  	} else {
  		interpolator = new irtkNearestNeighborInterpolateImageFunction;
  	}
  }
  interpolator->SetInput(&target);
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

  // Find out the maximum similarity frame.
  double currentvalue,maxvalue;
  int maxphase;

  maxvalue = 0;
  maxphase = 0;

  // Fill histogram
  for (t = 0; t < from.GetT(); t++){
	  for (z = 0; z < from.GetZ(); z++) {
		  for (y = 0; y < from.GetY(); y++) {
			  for (x = 0; x < from.GetX(); x++) {

				  val = from(x, y, z, t);

				  irtkPoint p(x, y, z);
				  // Transform point into world coordinates
				  from.ImageToWorld(p);
				  // Transform point
				  transformation->Transform(p);
				  // Transform point into image coordinates
				  target.WorldToImage(p);

				  // A bad thing might happen for the 2D case.
				  if ((target.GetZ() == 1) &&
					  (p._z > 0.5 || p._z < -0.5)){
						  cerr << "Transformed point outside plane of 2D target image." << endl;
						  exit(1);
				  }

				  // 2D and in plane but out of FoV.
				  if ((target.GetZ() == 1) &&
					  (p._x <= x1 || p._x >= x2 ||
					  p._y <= y1 || p._y >= y2))
					  continue;

				  // 3D and out of FoV.
				  if ((target.GetZ() > 1) &&
					  (p._x <= x1 || p._x >= x2 ||
					  p._y <= y1 || p._y >= y2 ||
					  p._z <= z1 || p._z >= z2))
					  continue;

				  // Should be able to interpolate if we've got this far.

				  val = interpolator->EvaluateInside(p._x, p._y, p._z);

				  histogram.AddSample(from(x, y, z, t), val);

			  }
		  }
	  }
	  currentvalue = histogram.NormalizedMutualInformation();
	  if(currentvalue > maxvalue){
		  maxvalue = currentvalue;
		  maxphase = t;
	  }
	  histogram.Reset();
  }

  //extract maxphase to output
  irtkRealImage output;
  irtkImageAttributes attr = from.GetImageAttributes();
  attr._t = 1; attr._dt = 1;
  output.Initialize(attr);

  for (z = 0; z < from.GetZ(); z++) {
	  for (y = 0; y < from.GetY(); y++) {
		  for (x = 0; x < from.GetX(); x++) {
			  output.PutAsDouble(x,y,z,from.GetAsDouble(x,y,z,maxphase));
		  }
	  }
  }

  output.Write(output_name);
  cout //<< //"closest temporal index is: " 
	  << maxphase << endl;
}
