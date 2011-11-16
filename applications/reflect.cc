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

#include <irtkFileToImage.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: reflect [in] [out]" << endl;
  cerr << "<-x/-y/-z>                                 \t Reflect the x/y/z-direction\n\n";
  cerr << "<-rx/-ry/-rz>                              \t Revert the x/y/z-direction\n\n";
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, x, y, z, revert_x, revert_y, revert_z;

  revert_x = false;
  revert_y = false;
  revert_z = false;
  x = false;
  y = false;
  z = false;

  if (argc < 4) {
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  irtkFileToImage *reader = irtkFileToImage::New(input_name);
  irtkBaseImage *image = reader->GetOutput();

  while (argc > 1) {
      ok = false;
      if (strcmp(argv[1], "-x") == 0) {
          x = true;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-y") == 0) {
          y = true;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-z") == 0) {
          z = true;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-rx") == 0) {
          revert_x = true;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-ry") == 0) {
          revert_y = true;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-rz") == 0) {
          revert_z = true;
          argc--;
          argv++;
          ok = true;
      } else if (!ok) {
          cerr << "Invalid option : " << argv[1] << endl;
          exit(1);
      }
  }

  if(x) image->ReflectX();
  if(y) image->ReflectY();
  if(z) image->ReflectZ();

  if (revert_x) {
      irtkImageAttributes tmpatr;
      tmpatr = image->GetImageAttributes();
      tmpatr._xaxis[0] = -tmpatr._xaxis[0];
      tmpatr._xaxis[1] = -tmpatr._xaxis[1];
      tmpatr._xaxis[2] = -tmpatr._xaxis[2];
      image->PutOrientation(tmpatr._xaxis,tmpatr._yaxis,tmpatr._zaxis);
  }

  if (revert_y) {
      irtkImageAttributes tmpatr;
      tmpatr = image->GetImageAttributes();
      tmpatr._yaxis[0] = -tmpatr._yaxis[0];
      tmpatr._yaxis[1] = -tmpatr._yaxis[1];
      tmpatr._yaxis[2] = -tmpatr._yaxis[2];
      image->PutOrientation(tmpatr._xaxis,tmpatr._yaxis,tmpatr._zaxis);
  }

  if (revert_z) {
      irtkImageAttributes tmpatr;
      tmpatr = image->GetImageAttributes();
      tmpatr._zaxis[0] = -tmpatr._zaxis[0];
      tmpatr._zaxis[1] = -tmpatr._zaxis[1];
      tmpatr._zaxis[2] = -tmpatr._zaxis[2];
      image->PutOrientation(tmpatr._xaxis,tmpatr._yaxis,tmpatr._zaxis);
  }


  // Write region
  image->Write(output_name);

  // Be nice
  delete image;
  delete reader;

  return 0;
}
