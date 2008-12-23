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

#include <irtkBiasField.h>

char *output;

void usage()
{
  cerr << "Usage: biasfield [image] [bias] [output] <options>" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double x, y, z;
  int i, j, k, ok, padding;

  if (argc < 4) {
    usage();
    exit(1);
  }

  // Input image
  irtkGreyImage image;
  image.Read(argv[1]);
  argc--;
  argv++;

  // Create bias field
  irtkBSplineBiasField *biasfield = new irtkBSplineBiasField;
  biasfield->Read(argv[1]);
  argc--;
  argv++;

  // Output file name
  output = argv[1];
  argc--;
  argv++;

  // Default parameters
  padding    = MIN_GREY;

  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  for (k = 0; k < image.GetZ(); k++) {
    for (j = 0; j < image.GetY(); j++) {
      for (i = 0; i < image.GetX(); i++) {
        if (image(i, j, k) != padding) {
          x = i;
          y = j;
          z = k;
          image.ImageToWorld(x, y, z);
          image(i, j, k) = round(biasfield->Bias(x, y, z));
        } else {
          image(i, j, k) = padding;
        }
      }
    }
  }

  image.Write(output);
}

