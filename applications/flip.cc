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

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: flip [in] [out] [-xy/-xz/-xt/-yz/-yz/-zt]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  irtkGreyImage image;

  if (argc != 4) {
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read input
  image.Read(input_name);

  if (strcmp("-xy", argv[1]) == 0) {
    image.FlipXY();
  } else {
    if (strcmp("-xz", argv[1]) == 0) {
      image.FlipXZ();
    } else {
      if (strcmp("-xt", argv[1]) == 0) {
        image.FlipXT();
      } else {
        if (strcmp("-yz", argv[1]) == 0) {
          image.FlipYZ();
        } else {
          if (strcmp("-yt", argv[1]) == 0) {
            image.FlipYT();
          } else {
            if (strcmp("-zt", argv[1]) == 0) {
              image.FlipZT();
            } else {
              usage();
            }
          }
        }
      }
    }
  }


  // Write region
  image.Write(output_name);

  return 0;
}
