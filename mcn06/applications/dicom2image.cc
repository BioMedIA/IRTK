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
  int ok, negate;
	double pos[3], origin[3];
  irtkImageAttributes attr;
  
  // Parse filenames
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Default
  negate = False;
  pos[0] = 0;
  pos[1] = 0;
  pos[2] = 0;
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = 0;
  
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dimensions") == 0)) {
      argc--;
      argv++;
      attr._x = atoi(argv[1]);
      argc--;
      argv++;
      attr._y = atoi(argv[1]);
      argc--;
      argv++;
      attr._z = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-size") == 0)) {
      argc--;
      argv++;
      attr._dx = atof(argv[1]);
      argc--;
      argv++;
      attr._dy = atof(argv[1]);
      argc--;
      argv++;
      attr._dz = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-orientation") == 0)) {
      argc--;
      argv++;
      attr._xaxis[0] = atof(argv[1]);
      argc--;
      argv++;
      attr._xaxis[1] = atof(argv[1]);
      argc--;
      argv++;
      attr._xaxis[2] = atof(argv[1]);
      argc--;
      argv++;
      attr._yaxis[0] = atof(argv[1]);
      argc--;
      argv++;
      attr._yaxis[1] = atof(argv[1]);
      argc--;
      argv++;
      attr._yaxis[2] = atof(argv[1]);
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
  	attr._xaxis[0] *= -1;
  	attr._xaxis[1] *= -1;
  	attr._xaxis[2] *= -1;
  	attr._yaxis[0] *= -1;
  	attr._yaxis[1] *= -1;
  	attr._yaxis[2] *= -1;
  }

	attr._zaxis[0] = attr._xaxis[1]*attr._yaxis[2] - attr._xaxis[2]*attr._yaxis[1];
	attr._zaxis[1] = attr._xaxis[2]*attr._yaxis[0] - attr._xaxis[0]*attr._yaxis[2];
	attr._zaxis[2] = attr._xaxis[0]*attr._yaxis[1] - attr._xaxis[1]*attr._yaxis[0];

  // Create image with specified voxel size and orientation
  irtkGreyImage image(attr);

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
  fread(image.GetPointerToVoxels(), sizeof(irtkGreyPixel), attr._x*attr._y*attr._z, fp);

  cout << "Writing to " << output_name << endl;
  image.Write(output_name);
}
