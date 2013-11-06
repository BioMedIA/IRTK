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
  cerr << "Usage: convert [in] [out] <options>\n\n";
  cerr << "where <options> can be one or more of the following:\n";
  cerr << "<-char|uchar|short|ushort|float|double>    \t Output voxel type\n";
  cerr << "<-minmax value value>                      \t Output min and max intensity\n";
  cerr << "<-saturate q0 q1>                      \t Saturate the image (typical usage: q0=0.01, q1=0.99)\n";  
  cerr << "Please note that IRTK will flip Analyze in the y-direction when the image \n";
  cerr << "is read and written (for historical reasons). This means that the coordinate \n";
  cerr << "system which IRTK uses for Analyze images is different from that used by other \n";
  cerr << "software Image Registration Toolkit (IRTK) such as SPM or FSL. Please use the NIFTI file format \n";
  cerr << "instead (preferred option) or use the reflect with -y flag before converting from or to\n";
  cerr << "Analyze file format.\n";
  exit(1);
}

int main(int argc, char **argv)
{
  double min, max;
  double q0, q1;
  int ok, minmax, saturate, image_type;

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

  // Parse remaining options
  minmax = false;
  min    = 0;
  max    = 0;
  image_type = IRTK_VOXEL_SHORT;

  while (argc > 1) {
      ok = false;
      if (strcmp(argv[1], "-char") == 0) {
          image_type = IRTK_VOXEL_CHAR;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-uchar") == 0) {
          image_type = IRTK_VOXEL_UNSIGNED_CHAR;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-short") == 0) {
          image_type = IRTK_VOXEL_SHORT;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-ushort") == 0) {
          image_type = IRTK_VOXEL_UNSIGNED_SHORT;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-float") == 0) {
          image_type = IRTK_VOXEL_FLOAT;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-double") == 0) {
          image_type = IRTK_VOXEL_DOUBLE;
          argc--;
          argv++;
          ok = true;
      } else if (strcmp(argv[1], "-minmax") == 0) {
          argc--;
          argv++;
          min = atof(argv[1]);
          argc--;
          argv++;
          max = atof(argv[1]);
          argc--;
          argv++;
          minmax = true;
          ok = true;
      } else if (strcmp(argv[1], "-saturate") == 0) {
          argc--;
          argv++;
          q0 = atof(argv[1]);
          argc--;
          argv++;
          q1 = atof(argv[1]);
          argc--;
          argv++;
          saturate = true;
          ok = true;          
      } else if (!ok) {
          cerr << "Invalid option : " << argv[1] << endl;
          exit(1);
      }
  }

  // Read image
  irtkGenericImage<irtkRealPixel> image(input_name);

  // Saturate
  if (saturate)
    image.Saturate(q0,q1);

  // Scale image
  if (minmax) {
    if (min >= max) {
      cerr << "Minimum value larger or equal to maximum value" << endl;
      exit(1);
    }
    image.PutMinMaxAsDouble(min, max);
  }

  // Convert image
  switch (image_type) {
    case IRTK_VOXEL_CHAR: {
        irtkGenericImage<char> output = image;
        output.Write(output_name);
      }
      break;
    case IRTK_VOXEL_UNSIGNED_CHAR: {
        irtkGenericImage<unsigned char> output = image;
        output.Write(output_name);
      }
      break;
    case IRTK_VOXEL_SHORT: {
        irtkGenericImage<short> output = image;
        output.Write(output_name);
      }
      break;
    case IRTK_VOXEL_UNSIGNED_SHORT: {
        irtkGenericImage<unsigned short> output = image;
        output.Write(output_name);
      }
      break;
    case IRTK_VOXEL_FLOAT: {
        irtkGenericImage<float> output = image;
        output.Write(output_name);
        break;
      }
    case IRTK_VOXEL_DOUBLE: {
        irtkGenericImage<double> output = image;
        output.Write(output_name);
        break;
      }
    default:
      cerr << "Unknown voxel type for output format" << endl;
      exit(1);
  }
}
