/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: headertool.cc 235 2010-10-18 09:25:20Z dr $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2010-10-18 10:25:20 +0100 (一, 18 十月 2010) $
  Version   : $Revision: 235 $
  Changes   : $Author: dr $

=========================================================================*/

#include <irtkImage.h>

#include <irtkFileToImage.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: headertool [in] [out] <options>\n" << endl;
  cerr << "where <options> can be one or more of the following:\n";
  cerr << "<-size        dx dy dz>                      \t Voxel size   (in mm)\n";
  cerr << "<-tsize       dt>                            \t Voxel size   (in ms)\n";
  cerr << "<-orientation x1 x2 x3   y1 y2 y3  z1 z2 z3> \t Image orientation\n";
  cerr << "<-origin      x  y  z>                       \t Image origin (in mm)\n";
  cerr << "<-torigin     t>                             \t Image origin (in ms)\n";
  cerr << "<-dofin       dof>                           \t Change header from transformation\n";
  cerr << "<-target      image>                         \t Copy header from target image\n\n";
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  double xsize, ysize, zsize, tsize, xaxis[3], yaxis[3], zaxis[3], origin[4];

  if (argc < 3) {
    usage();
  }

  // Parse filenames
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read image
  irtkFileToImage *reader = irtkFileToImage::New(input_name);
  irtkBaseImage *image = reader->GetOutput();

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-target") == 0)) {
      argc--;
      argv++;
      irtkGreyImage target(argv[1]);
      argc--;
      argv++;
      target.GetPixelSize(&xsize, &ysize, &zsize);
      image->PutPixelSize(xsize, ysize, zsize);
      target.GetOrientation(xaxis, yaxis, zaxis);
      image->PutOrientation(xaxis, yaxis, zaxis);
      target.GetOrigin(origin[0], origin[1], origin[2]);
      image->PutOrigin(origin[0], origin[1], origin[2]);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-size") == 0)) {
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
      ok = true;

      // Set voxel size
      image->PutPixelSize(xsize, ysize, zsize);
    }
    if ((ok == false) && (strcmp(argv[1], "-tsize") == 0)) {
      argc--;
      argv++;
      tsize = atof(argv[1]);
      argc--;
      argv++;
      ok = true;

      // Set voxel size
      image->GetPixelSize(&xsize, &ysize, &zsize);
      image->PutPixelSize(xsize, ysize, zsize, tsize);
    }
    if ((ok == false) && (strcmp(argv[1], "-orientation") == 0)) {
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
      zaxis[0] = atof(argv[1]);
      argc--;
      argv++;
      zaxis[1] = atof(argv[1]);
      argc--;
      argv++;
      zaxis[2] = atof(argv[1]);
      argc--;
      argv++;
      ok = true;

      // Set orientation (third argument now required)
      image->PutOrientation(xaxis, yaxis, zaxis);
    }
    if ((ok == false) && (strcmp(argv[1], "-origin") == 0)) {
      argc--;
      argv++;
      origin[0] = atof(argv[1]);
      argc--;
      argv++;
      origin[1] = atof(argv[1]);
      argc--;
      argv++;
      origin[2] = atof(argv[1]);
      argc--;
      argv++;
      ok = true;

      // Set origin
      image->PutOrigin(origin[0], origin[1], origin[2]);
    }
    if ((ok == false) && (strcmp(argv[1], "-torigin") == 0)) {
      argc--;
      argv++;
      origin[3] = atof(argv[1]);
      argc--;
      argv++;
      ok = true;

      // Set origin
      image->GetOrigin(origin[0], origin[1], origin[2]);
      image->PutOrigin(origin[0], origin[1], origin[2], origin[3]);
    }

    if (ok == false) {
      cout << "Can't parse argument: " << argv[1] << endl;
      usage();
    }
  }

  image->Write(output_name);

  return 0;
}

