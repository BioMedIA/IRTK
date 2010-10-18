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

#include <irtkGaussianBlurring.h>

#include <irtkGaussianBlurring4D.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: blur [in] [out] [sigma] <options>" << endl;
  cerr << "where <options> are one or more of the following:\n";
  cerr << "\t<-3D>              3D blurring (default)\n";
  cerr << "\t<-4D>              4D blurring\n";
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, blur4D;
  double sigma;
  irtkGreyImage input;

  if (argc < 4) {
    usage();
  }

  // Parse parameters
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  sigma = atof(argv[1]);
  argc--;
  argv++;

  // Default
  blur4D = false;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-3D") == 0)) {
      argc--;
      argv++;
      blur4D = false;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-4D") == 0)) {
      argc--;
      argv++;
      blur4D = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  // Read input
  input.Read(input_name);

  // Blur image
  if (blur4D == true) {
    irtkGaussianBlurring4D<irtkGreyPixel> gaussianBlurring4D(sigma);
    gaussianBlurring4D.SetInput (&input);
    gaussianBlurring4D.SetOutput(&input);
    gaussianBlurring4D.Run();
  } else {
    irtkGaussianBlurring<irtkGreyPixel> gaussianBlurring(sigma);
    gaussianBlurring.SetInput (&input);
    gaussianBlurring.SetOutput(&input);
    gaussianBlurring.Run();
  }

  // Write image
  input.Write(output_name);

  return 0;
}
