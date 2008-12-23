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
#include <irtkResampling.h>
#include <irtkResamplingWithPadding.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: resample [in] [out] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "\t<-size x y z>      New voxel size in mm (default: 1 1 1)\n";
  cerr << "\t<-linear>          Linear interpolation\n";
  cerr << "\t<-bspline>         B-spline interpolation\n";
  cerr << "\t<-cspline>         Cubic spline interpolation\n";
  cerr << "\t<-sinc>            Truncated sinc interpolation\n";
  cerr << "\t<-gaussian sigma>  Gaussian interpolation\n";
  cerr << "\t<-padding value>   Background padding (default: MIN_SHRT)\n";
  cerr << "\t                   Only linear interpolation for padding!\n\n";
  exit(1);
}

int main(int argc, char **argv)
{
  Bool ok, padding;
  double xsize, ysize, zsize;
  irtkGreyImage image;
  irtkImageFunction<irtkGreyPixel> *interpolator = NULL;
  irtkGreyPixel padding_value = MIN_GREY;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Read input and output names
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read input
  image.Read(input_name);

  // Parse remaining parameters
  xsize = 1;
  ysize = 1;
  zsize = 1;
  padding = False;
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-size") == 0)) {
      argc--;
      argv++;
      xsize = atof(argv[1]);
      argc--;
      argv++;
      ysize = atof(argv[1]);
      argc--;
      argv++;
      zsize = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-linear") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkLinearInterpolateImageFunction<irtkGreyPixel>;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-bspline") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkBSplineInterpolateImageFunction<irtkGreyPixel>;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-cspline") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkCSplineInterpolateImageFunction<irtkGreyPixel>;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-sinc") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkSincInterpolateImageFunction<irtkGreyPixel>;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-gaussian") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkGaussianInterpolateImageFunction<irtkGreyPixel>(atof(argv[1]));
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding_value = atoi(argv[1]);
      argc--;
      argv++;
      padding = True;
      ok = True;
    }
    if (ok == False) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  // Create default interpolator
  if (interpolator == NULL) {
    interpolator = new irtkNearestNeighborInterpolateImageFunction<irtkGreyPixel>;
  }

  // Resample
  if (padding == False) {
    cout << "Resampling ... "; cout.flush();
    irtkResampling<irtkGreyPixel> resampling(xsize, ysize, zsize);
    resampling.SetInput(&image);
    resampling.SetOutput(&image);
    resampling.SetInterpolator(interpolator);
    resampling.Run();
    cout << "done" << endl;
  } else {
    cout << "Resampling with padding ... "; cout.flush();
    irtkResamplingWithPadding<irtkGreyPixel> resampling(xsize, ysize, zsize,
        padding_value);
    resampling.SetInput(&image);
    resampling.SetOutput(&image);
    resampling.SetInterpolator(interpolator);
    resampling.Run();
    cout << "done" << endl;
  }

  // Save result
  image.Write(output_name);

  return 0;
}
