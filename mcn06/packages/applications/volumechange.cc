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
char *input_name = NULL, *output_name, *dof_name  = NULL;

void usage()
{
  cerr << "Usage: volumechange [mask] [mask-value] [dof]\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, n, mask_value;
  double x, y, z, volume;

  // Check command line
  if (argc != 4) {
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  mask_value = atoi(argv[1]);
  argc--;
  argv++;
  dof_name = argv[1];
  argc--;
  argv++;

  // Read image
  cout << "Reading image ... "; cout.flush();
  irtkGreyImage *image = new irtkGreyImage(input_name);
  cout << "done" << endl;

  // Read transformation
  irtkTransformation *transformation = irtkTransformation::New(dof_name);

  n      = 0;
  volume = 0;
  for (k = 0; k < image->GetZ(); k++) {
    for (j = 0; j < image->GetY(); j++) {
      for (i = 0; i < image->GetX(); i++) {
        if (image->Get(i, j, k) == mask_value) {
          x = i;
          y = j;
          z = k;
          image->ImageToWorld(x, y, z);
          volume += transformation->Jacobian(x, y, z);
          n++;
        }
      }
    }
  }
  cout << "Volume change is " << ((volume / n) - 1) * 100 << "%" << endl;
}
