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
  cerr << "Usage: binarize [in] [out] [threshold] [value1] [value2]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i;
  irtkGreyImage input;
  irtkGreyPixel *ptr, threshold, value1, value2;

  if (argc != 6) {
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  threshold = atoi(argv[1]);
  argc--;
  argv++;
  value1 = atoi(argv[1]);
  argc--;
  argv++;
  value2 = atoi(argv[1]);
  argc--;
  argv++;

  // Read input
  input.Read(input_name);

  ptr = input.GetPointerToVoxels();
  for (i = 0; i < input.GetNumberOfVoxels(); i++) {
    if (*ptr <= threshold) {
      *ptr = value1;
    } else {
      *ptr = value2;
    }
    ptr++;
  }

  // Write image
  input.Write(output_name);

  return 0;
}
