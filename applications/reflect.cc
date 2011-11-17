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
  cerr << " Usage: reflect [in] [out] [reflection_1] <reflection_2> .. <reflection_n>" << endl;
  cerr << " " << endl;
  cerr << " Apply a sequence of one or more spatial reflections to an image." << endl;
  cerr << " Where the reflections are chosen from:" << endl;
  cerr << " " << endl;
  cerr << "        -x -y -z -xy -xz -yz" << endl;
  cerr << " " << endl;
  cerr << " The reflections are processed in the order given" << endl;
  cerr << " and changing the order can change the result. " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;

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

      if ((ok == false) && (strcmp(argv[1], "-x") == 0)) {
          argc--;
          argv++;
          image->ReflectX();
          ok = true;
      }
      if ((ok == false) && (strcmp(argv[1], "-y") == 0)) {
          argc--;
          argv++;
          image->ReflectY();
          ok = true;
      }
      if ((ok == false) && (strcmp(argv[1], "-z") == 0)) {
          argc--;
          argv++;
          image->ReflectZ();
          ok = true;
      }
      if ((ok == false) && ((strcmp(argv[1], "-xy") == 0) || (strcmp(argv[1], "-yx") == 0))) {
          argc--;
          argv++;
          image->FlipXY(0);
          ok = true;
      }
      if ((ok == false) && ((strcmp(argv[1], "-xz") == 0) || (strcmp(argv[1], "-zx") == 0))) {
          argc--;
          argv++;
          image->FlipXZ(0);
          ok = true;
      }
      if ((ok == false) && ((strcmp(argv[1], "-yz") == 0) || (strcmp(argv[1], "-zy") == 0))) {
          argc--;
          argv++;
          image->FlipYZ(0);
          ok = true;
      }

      if (ok == false) {
  			cout << "Can't parse argument: " << argv[1] << endl;
  			usage();
      }
  }

  // Write image
  image->Write(output_name);

  // Be nice
  delete image;
  delete reader;

  return 0;
}
