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

#include <irtkGradientImage.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: edgedetect [in] [out] [sigma]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkRealImage input;

  if (argc < 4) {
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read input
  input.Read(input_name);

  // Blur image
  if(atof(argv[1]) > 0){
      irtkGaussianBlurring<irtkRealPixel> gaussianBlurring(atof(argv[1]));
      gaussianBlurring.SetInput (&input);
      gaussianBlurring.SetOutput(&input);
      gaussianBlurring.Run();
  }

  // Compute gradient
  irtkGradientImage<irtkRealPixel> gradient;
  gradient.SetInput (&input);
  gradient.SetOutput(&input);
  gradient.Run();

  // Write image
  input.Write(output_name);

  return 0;
}
