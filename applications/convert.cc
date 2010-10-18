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
  cerr << "where <options can be one or more of the following:\n";
  cerr << "<-char|uchar|short|ushort|float|double>    Output voxel type\n";
  cerr << "<-minmax value value>                      Output min and max intensity\n";
  cerr << "<-x/-y/-z>                                 Flip the image in the x/y/z-direction\n\n";
  cerr << "Please note that IRTK will flip Analyze in the y-direction when the image \n";
  cerr << "is read and written (for historical reasons). This means that the coordinate \n";
  cerr << "system which IRTK uses for Analyze images is different from that used by other \n";
  cerr << "software Image Registration Toolkit (IRTK) such as SPM or FSL. Please use the NIFTI file format \n";
  cerr << "instead (preferred option) or use the -y flag before converting from or to\n";
  cerr << "Analyze file format.\n";
  exit(1);
}

int main(int argc, char **argv)
{
  double min, max;
  int ok, minmax, flip_x, flip_y, flip_z, image_type;

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
  flip_x = false;
  flip_y = false;
  flip_z = false;
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
    } else if (strcmp(argv[1], "-x") == 0) {
      flip_x = true;
      argc--;
      argv++;
      ok = true;
    } else if (strcmp(argv[1], "-y") == 0) {
      flip_y = true;
      argc--;
      argv++;
      ok = true;
    } else if (strcmp(argv[1], "-z") == 0) {
      flip_y = true;
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
    } else if (!ok) {
      cerr << "Invalid option : " << argv[1] << endl;
      exit(1);
    }
  }
  
  // Read image
  irtkGenericImage<double> image(input_name);

  // Scale image
  if (minmax) {
    if (min >= max) {
      cerr << "Minimum value larger or equal to maximum value" << endl;
      exit(1);
    }
    image.PutMinMaxAsDouble(min, max);
  }
  
  // Reflect image
  if (flip_x == true) image.ReflectX();
  if (flip_y == true) image.ReflectY();
  if (flip_z == true) image.ReflectZ();

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
