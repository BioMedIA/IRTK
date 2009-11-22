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
  cerr << "Usage: reflect [in] [out] [-x/-y/-z]" << endl;
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

  if (strcmp("-x", argv[1]) == 0) {
    image.ReflectX();
  } else {
    if (strcmp("-y", argv[1]) == 0) {
      image.ReflectY();
    } else {
      if (strcmp("-z", argv[1]) == 0) {
        image.ReflectZ();
      } else {
        usage();
      }
    }
  }


  // Write region
  image.Write(output_name);

  return 0;
}
