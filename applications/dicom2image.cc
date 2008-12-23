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
  cerr << "Usage: dicom2image [in] [out] <options>\n" << endl;
  exit(1);

}

int main(int argc, char **argv)
{
  int ok, negate, x, y, z;
  double xsize, ysize, zsize, xaxis[3], yaxis[3], origin[3], pos[3];

  // Parse filenames
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Default image size
  x = 1;
  y = 1;
  z = 1;

  // Default voxel size
  xsize = 1;
  ysize = 1;
  zsize = 1;

  // Default orientation
  xaxis[0] = 1;
  xaxis[1] = 0;
  xaxis[2] = 0;
  yaxis[0] = 0;
  yaxis[1] = 1;
  yaxis[2] = 0;

  // Default position
  pos[0] = 0;
  pos[1] = 0;
  pos[2] = 0;

  // Default
  negate = False;

  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dimensions") == 0)) {
      argc--;
      argv++;
      x = atoi(argv[1]);
      argc--;
      argv++;
      y = atoi(argv[1]);
      argc--;
      argv++;
      z = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-size") == 0)) {
      argc--;
      argv++;
      xsize = atof(argv[1]);
      argc--;
      argv++;
      ysize = atof(argv[1]);
      argc--;
      argv++;
      zsize = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-orientation") == 0)) {
      argc--;
      argv++;
      xaxis[0] = atof(argv[1]);
      argc--;
      argv++;
      xaxis[1] = atof(argv[1]);
      argc--;
      argv++;
      xaxis[2] = atof(argv[1]);
      argc--;
      argv++;
      yaxis[0] = atof(argv[1]);
      argc--;
      argv++;
      yaxis[1] = atof(argv[1]);
      argc--;
      argv++;
      yaxis[2] = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-position") == 0)) {
      argc--;
      argv++;
      pos[0] = atof(argv[1]);
      argc--;
      argv++;
      pos[1] = atof(argv[1]);
      argc--;
      argv++;
      pos[2] = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-negate") == 0)) {
      argc--;
      argv++;
      negate = True;
      ok = True;
    }
    if (ok == False) {
      cout << "Can't parse argument: " << argv[1] << endl;
      usage();
    }
  }

  // Negate coordinate axis
  if (negate == True) {
    xaxis[0] *= -1;
    xaxis[1] *= -1;
    xaxis[2] *= -1;
    yaxis[0] *= -1;
    yaxis[1] *= -1;
    yaxis[2] *= -1;
  }

  // Create i\mage with specified voxel size and orientation
  irtkGreyImage image(x, y, z, xsize, ysize, zsize, xaxis, yaxis);

  // Calculate position of first voxel
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = 0;
  image.ImageToWorld(origin[0], origin[1], origin[2]);

  // Calculate new origin
  origin[0] = pos[0] - origin[0];
  origin[1] = pos[1] - origin[1];
  origin[2] = pos[2] - origin[2];
  image.PutOrigin(origin[0], origin[1], origin[2]);

  // Read raw data
  FILE *fp = fopen(input_name, "r");
  fread(image.GetPointerToVoxels(), sizeof(irtkGreyPixel), x*y*z, fp);

  cout << "Writing to " << output_name << endl;
  image.Write(output_name);
}
