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
  cerr << "Usage: padding [imageA] [imageB] [output] [Value in imageB] [Padding Value in output]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z;
  irtkGreyPixel threshold, padding;

  if (argc != 6) {
    usage();
  }

  irtkGreyImage imageA(argv[1]);
  irtkGreyImage imageB(argv[2]);

  threshold = atoi(argv[4]);
  padding = atoi(argv[5]);

  for (z = 0; z < imageA.GetZ(); z++) {
    for (y = 0; y < imageA.GetY(); y++) {
      for (x = 0; x < imageA.GetX(); x++) {
        if (imageB(x, y, z) == threshold) imageA(x, y, z) = padding;
      }
    }
  }
  imageA.Write(argv[3]);

  return 0;
}
