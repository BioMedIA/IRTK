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

int main(int argc, char **argv)
{
  irtkGreyImage image;
  double xsize, ysize, zsize;

  if (argc != 6) {
    cerr << "usage: voxelsize [input] [output] [xsize] [ysize] [zsize]" << endl;
    exit(1);
  }
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;
  xsize = atof(argv[1]);
  argc--;
  argv++;
  ysize = atof(argv[1]);
  argc--;
  argv++;
  zsize = atof(argv[1]);

  image.Read(input_name);
  image.PutPixelSize(xsize, ysize, zsize);
  image.Write(output_name);

  return 0;
}

