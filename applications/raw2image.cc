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

void usage()
{
  cerr << "Usage: raw2image [rawfile] [image] [x] [y] [z] [float|short|ushort|char|uchar]\n";
  exit(1);
}

int main(int argc, char **argv)
{
  FILE *fp;
  char *input_name, *output_name;
  int x, y, z;

  if (argc < 6) {
    usage();
  }

  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;
  x = atoi(argv[1]);
  argv++;
  argc--;
  y = atoi(argv[1]);
  argv++;
  argc--;
  z = atoi(argv[1]);
  argv++;
  argc--;

  size_t nvoxels = x*y*z;

  // Open raw file
  fp = fopen(input_name, "r");
  if (fp == NULL) {
    perror("fopen");
    exit(1);
  }

  if (strcmp(argv[1], "float") == 0) {
    irtkRealImage image(x, y, z);
    if (fread(image.GetPointerToVoxels(), sizeof(float), nvoxels, fp) != nvoxels) {
      fclose(fp);
      cerr << "Failed to read data from input file or file contains less voxels than specified" << endl;
      exit(1);
    }
    image.Write(output_name);
  } else {
    if (strcmp(argv[1], "short") == 0) {
      irtkGreyImage image(x, y, z);
      if (fread(image.GetPointerToVoxels(), sizeof(short), nvoxels, fp) != nvoxels) {
        fclose(fp);
        cerr << "Failed to read data from input file or file contains less voxels than specified" << endl;
        exit(1);
      }
      image.Write(output_name);
    } else {
      if (strcmp(argv[1], "uchar") == 0) {
        irtkByteImage image(x, y, z);
        if (fread(image.GetPointerToVoxels(), sizeof(unsigned char), nvoxels, fp) != nvoxels) {
          fclose(fp);
          cerr << "Failed to read data from input file or file contains less voxels than specified" << endl;
          exit(1);
        }
        image.Write(output_name);
      } else {
        cerr << "Format: " << argv[1] << " not yet supported" << endl;
        exit(1);
      }
    }
  }

  fclose(fp);
  return 0;
}
