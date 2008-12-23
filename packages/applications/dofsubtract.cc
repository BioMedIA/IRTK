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

#include <irtkTransformation.h>

// Default filenames
char *image_name = NULL, *input1_name = NULL, *input2_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: dof2subtract [image] [dof1] [dof2] [output] ";
  cerr << "<options> \n\nwhere <options> is one or more of the ";
  cerr << "following:\n\n";
  cerr << "-padding value\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z, ok, padding;
  double p1[3], p2[3];

  if (argc < 4) {
    usage();
  }

  // Parse file names
  image_name  = argv[1];
  argc--;
  argv++;
  input1_name  = argv[1];
  argc--;
  argv++;
  input2_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Parse remaining parameters
  padding = MIN_GREY;
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

  // Read image
  irtkGreyImage image(image_name);

  // Read transformation
  irtkTransformation *transform1 = irtkTransformation::New(input1_name);

  // Read transformation
  irtkTransformation *transform2 = irtkTransformation::New(input2_name);

  // Initialize point structure with transformed point positions
  for (z = 0; z < image.GetZ(); z++) {
    for (y = 0; y < image.GetY(); y++) {
      for (x = 0; x < image.GetX(); x++) {
        if (image(x, y, z) > padding) {
          p1[0] = x;
          p1[1] = y;
          p1[2] = z;
          image.ImageToWorld(p1[0], p1[1], p1[2]);
          p2[0] = p1[0];
          p2[1] = p1[1];
          p2[2] = p1[2];
          transform1->Transform(p1[0], p1[1], p1[2]);
          transform2->Transform(p2[0], p2[1], p2[2]);
          image(x, y, z) = round(100 * sqrt(pow(p1[0] - p2[0], 2.0) +
                                            pow(p1[1] - p2[1], 2.0) +
                                            pow(p1[2] - p2[2], 2.0)));
        } else {
          image(x, y, z) = 0;
        }
      }
    }
  }
  image.Write(output_name);

}

