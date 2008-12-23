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

typedef enum pixel_type {pt_real, pt_grey, pt_byte, pt_rgb} pixel_type;

void usage()
{
  cerr << "Usage: convert [in] [out] <options>\n\n";
  cerr << "\t where <options can be one or more of the following:\n";
  cerr << "\t\t <-real|-grey|-byte|-rgb>    Output voxel type\n";
  cerr << "\t\t <-minmax value value>       Output min and max intensity\n";
  cerr << "\t\t <-x/-y/-z>                  Flip the image in the x/y/z-direction\n\n";
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
  Bool ok, minmax, flip_x, flip_y, flip_z;
  pixel_type pt = pt_grey;
  float min = 0, max = 0;

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
  flip_x = False;
  flip_y = False;
  flip_z = False;
  minmax = False;
  while (argc > 1) {
    ok = False;
    if (strcmp(argv[1], "-real") == 0) {
      pt = pt_real;
      argc--;
      argv++;
      ok = True;
    } else if (strcmp(argv[1], "-grey") == 0) {
      pt = pt_grey;
      argc--;
      argv++;
      ok = True;
    } else if (strcmp(argv[1], "-byte") == 0) {
      pt = pt_byte;
      argc--;
      argv++;
      ok = True;
    } else if (strcmp(argv[1], "-rgb") == 0) {
      pt = pt_rgb;
      argc--;
      argv++;
      ok = True;
    } else if (strcmp(argv[1], "-x") == 0) {
      flip_x = True;
      argc--;
      argv++;
      ok = True;
    } else if (strcmp(argv[1], "-y") == 0) {
      flip_y = True;
      argc--;
      argv++;
      ok = True;
    } else if (strcmp(argv[1], "-z") == 0) {
      flip_y = True;
      argc--;
      argv++;
      ok = True;
    } else if (strcmp(argv[1], "-minmax") == 0) {
      argc--;
      argv++;
      min = atof(argv[1]);
      argc--;
      argv++;
      max = atof(argv[1]);
      argc--;
      argv++;
      minmax = True;
      ok = True;
    } else if (!ok) {
      cerr << "Invalid option : " << argv[1] << endl;
      exit(1);
    }
  }

  // Save output
  switch (pt) {
  case pt_real: {
    // Read in as irtkRealImage
    irtkRealImage image(input_name);

    if (minmax) {
      if (min >= max) {
        cerr << "Minimum value larger or equal to maximum value" << endl;
        exit(1);
      }
      image.PutMinMax(min, max);
    }
    if (flip_x == True) image.ReflectX();
    if (flip_y == True) image.ReflectY();
    if (flip_z == True) image.ReflectZ();
    image.Write(output_name);
  }
  break;
  case pt_grey: {
    // Read in as irtkGreyImage
    irtkGreyImage image(input_name);

    if (minmax) {
      irtkGreyPixel min_val = round(min);
      irtkGreyPixel max_val = round(max);
      if (min < MIN_GREY) min_val = 0;
      if (max > MAX_GREY) max_val = MAX_GREY;
      if (min_val >= max_val) {
        cerr << "Minimum value larger or equal to maximum value" << endl;
        exit(1);
      }
      image.PutMinMax(min_val, max_val);
    }
    if (flip_x == True) image.ReflectX();
    if (flip_y == True) image.ReflectY();
    if (flip_z == True) image.ReflectZ();
    image.Write(output_name);
  }
  break;
  case pt_byte: {
    // Read in as irtkByteImage
    irtkByteImage image(input_name);

    if (minmax) {
      irtkBytePixel min_val = round(min);
      irtkBytePixel max_val = round(max);
      if (min < 0) min_val = 0;
      if (max > MAX_BYTE) max_val = MAX_BYTE;
      if (min_val >= max_val) {
        cerr << "Minimum value larger or equal to maximum value" << endl;
      }
      image.PutMinMax(min_val, max_val);
    }
    if (flip_x == True) image.ReflectX();
    if (flip_y == True) image.ReflectY();
    if (flip_z == True) image.ReflectZ();
    image.Write(output_name);
  }
  break;
  case pt_rgb: {
    // Read in as irtkRealImage
    irtkRealImage image(input_name);

    // Convert to irtkRGBImage
    irtkRGBImage image2 = image;

    // Write image
    if (flip_x == True) image.ReflectX();
    if (flip_y == True) image.ReflectY();
    if (flip_z == True) image.ReflectZ();
    image2.Write(output_name);
  }
  break;
  }

  return 0;
}
