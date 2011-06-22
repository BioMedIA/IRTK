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

typedef enum pixel_type {pt_real, pt_grey, pt_byte} pixel_type;

void usage()
{
    cerr << "Usage: rescale [in] [out] <min> <max>\n";
    exit(1);
}

int main(int argc, char **argv)
{
  irtkGreyPixel min, max;

  if (argc != 5) {
    usage();
  }

  // Parse filenames
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  min = atoi(argv[1]);
  argc--;
  argv++;
  max = atoi(argv[1]);
  argc--;
  argv++;

  // Read
  irtkGreyImage image(input_name);
  image.PutMinMax(min, max);
  image.Write(output_name);

  return 0;
}
