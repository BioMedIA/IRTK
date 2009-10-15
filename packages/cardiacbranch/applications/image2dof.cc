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
char *input_name = NULL, *input_name_x = NULL, *input_name_y = NULL, *input_name_z = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: image2dof [image] [dx] [dy] [dz] [dof]\n" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int x, y, z;
  double x1, y1, z1, x2, y2, z2, xsize, ysize, zsize, xaxis[3], yaxis[3], zaxis[3];

  if (argc < 3) {
    usage();
  }

  // Parse file names
  input_name   = argv[1];
  argc--;
  argv++;
  input_name_x = argv[1];
  argc--;
  argv++;
  input_name_y = argv[1];
  argc--;
  argv++;
  input_name_z = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;

  // Read image
  irtkGreyImage image;
  image.Read(input_name);

  // Read images with displacements
  irtkRealImage dx;
  irtkRealImage dy;
  irtkRealImage dz;
  dx.Read(input_name_x);
  dy.Read(input_name_y);
  dz.Read(input_name_z);

  // Compute bounding box of data
  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = dx.GetX()-1;
  y2 = dx.GetY()-1;
  z2 = dx.GetZ()-1;
  dx.ImageToWorld(x1, y1, z1);
  dx.ImageToWorld(x2, y2, z2);
  dx.GetOrientation(xaxis, yaxis, zaxis);
  dx.GetPixelSize(&xsize, &ysize, &zsize);

  // Create transformation
  irtkMultiLevelFreeFormTransformation *transformation = new irtkMultiLevelFreeFormTransformation;

  // Create deformation
  irtkLinearFreeFormTransformation *ffd = new
  irtkLinearFreeFormTransformation(x1, y1, z1, x2, y2, z2, xsize, ysize, zsize, xaxis, yaxis, zaxis);

  // Initialize point structure with transformed point positions
  for (z = 0; z < dx.GetZ(); z++) {
    for (y = 0; y < dx.GetY(); y++) {
      for (x = 0; x < dx.GetX(); x++) {
        ffd->Put(x, y, z, dx(x, y, z), dy(x, y, z), dz(x, y, z));
      }
    }
  }

  // Add deformation
  transformation->PushLocalTransformation(ffd);

  // Write transformation
  transformation->irtkTransformation::Write(output_name);
}
