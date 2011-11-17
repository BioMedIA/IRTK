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
  cerr << " Usage: flip [in] [out] [-xy/-xz/-xt/-yz/-yz/-zt]" << endl;
  cerr << " " << endl;
  cerr << " Swap the specified dimensions of the image grid _and_ " << endl;
  cerr << " swap the corresponding world coordinates of the centre " << endl;
  cerr << " of the grid (origin information in header)." << endl;
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
    image.FlipXY(1);
  } else {
    if (strcmp("-xz", argv[1]) == 0) {
      image.FlipXZ(1);
    } else {
      if (strcmp("-xt", argv[1]) == 0) {
        image.FlipXT(1);
      } else {
        if (strcmp("-yz", argv[1]) == 0) {
          image.FlipYZ(1);
        } else {
          if (strcmp("-yt", argv[1]) == 0) {
            image.FlipYT(1);
          } else {
            if (strcmp("-zt", argv[1]) == 0) {
              image.FlipZT(1);
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
