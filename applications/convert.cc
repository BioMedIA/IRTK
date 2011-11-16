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
  cerr << "<-swapzt>								  \t Swap z t axis but preserve the origin of z (for xy and so on use flip)\n";
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
  int ok, i,j,k,t, minmax, image_type, swapzt;

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
  swapzt = 0;

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
      } else if (strcmp(argv[1], "-swapzt") == 0) {
          argc--;
          argv++;
          swapzt = 1;
          ok = true;
      } else if (!ok) {
          cerr << "Invalid option : " << argv[1] << endl;
          exit(1);
      }
  }

  // Read image
  irtkGenericImage<float> image(input_name);

  // Scale image
  if (minmax) {
    if (min >= max) {
      cerr << "Minimum value larger or equal to maximum value" << endl;
      exit(1);
    }
    image.PutMinMaxAsDouble(min, max);
  }

  if(swapzt == 1) {
    irtkImageAttributes atrx;
    atrx = image.GetImageAttributes();
    irtkImageAttributes atry;
    atry = image.GetImageAttributes();
    atry._z = atrx._t;
    atry._t = atrx._z;
    atry._dz = atrx._dz;
    atry._dt = atrx._dt;
    if (atry._dt != 1) {
      atry._dt = 1;
    }
    image.WorldToImage(atry._xorigin,atry._yorigin,atry._zorigin);
    atry._zorigin -= atrx._z / 2.0;
    atry._zorigin += 0.5;
    image.ImageToWorld(atry._xorigin,atry._yorigin,atry._zorigin);

    irtkGreyImage toutput(atry);
    for(t=0;t<atrx._t;t++) {
      for(k = 0; k < atrx._z; k++) {
        for(j = 0; j < atrx._y; j++) {
          for(i = 0; i < atrx._x; i++) {
            toutput.PutAsDouble(i,j,t,k,image.GetAsDouble(i,j,k,t));
          }
        }
      }
    }
    image = toutput;
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
