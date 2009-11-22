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

void usage()
{
  cerr << "Usage: threshold [input] [output] [threshold]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z;
  irtkGreyPixel threshold;

  if (argc != 4) {
    usage();
  }

  // Read image
  irtkGreyImage image(argv[1]);

  // Threshold value
  threshold = atoi(argv[3]);

  for (z = 0; z < image.GetZ(); z++) {
    for (y = 0; y < image.GetY(); y++) {
      for (x = 0; x < image.GetX(); x++) {
        if (image(x, y, z) <= threshold) {
          image(x, y, z) = 0;
        } else {
          image(x, y, z) = 1;
        }
      }
    }
  }

  image.Write(argv[2]);

  return 0;
}
