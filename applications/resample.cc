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
#include <irtkImageFunction.h>
#include <irtkResampling.h>
#include <irtkResamplingWithPadding.h>

char *input_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: resample [in] [out] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "\t<-size x y z>      New voxel size in mm (default: 1 1 1)\n";
  cerr << "\t<-isotropic m>     Resample to isotropic size z size*m\n";
  cerr << "\t<-linear>          Linear interpolation\n";
  cerr << "\t<-bspline>         B-spline interpolation\n";
  cerr << "\t<-cspline>         Cubic spline interpolation\n";
  cerr << "\t<-sinc>            Truncated sinc interpolation\n";
  cerr << "\t<-sbased>          Shape based interpolation\n";
  cerr << "\t<-gaussian sigma>  Gaussian interpolation\n";
  cerr << "\t<-padding value>   Background padding (default: MIN_SHRT)\n";
  cerr << "\t                   Only linear interpolation for padding!\n\n";
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok, padding;
  double xsize, ysize, zsize, isotropic;
  irtkImage *image;
  irtkImageFunction *interpolator = NULL;
  irtkGreyPixel padding_value = MIN_GREY;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Read input and output names
  input_name  = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Read input
  image = irtkImage::New(input_name);

  // Parse remaining parameters
  xsize = 1;
  ysize = 1;
  zsize = 1;
  padding = false;
  isotropic = 0;

  while (argc > 1) {
    ok = false;
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
    }
    if ((ok == false) && (strcmp(argv[1], "-linear") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkLinearInterpolateImageFunction;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-isotropic") == 0)) {
      argc--;
      argv++;
      if (argc > 1 && argv[1][0] != '-') {
        isotropic = atof(argv[1]);
        argc--;
        argv++;
      } else
        isotropic = 1.0;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkBSplineInterpolateImageFunction;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-cspline") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkCSplineInterpolateImageFunction;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sinc") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkSincInterpolateImageFunction;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sbased") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkShapeBasedInterpolateImageFunction;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-gaussian") == 0)) {
      argc--;
      argv++;
      interpolator = new irtkGaussianInterpolateImageFunction(atof(argv[1]));
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-padding") == 0)) {
      argc--;
      argv++;
      padding_value = atoi(argv[1]);
      argc--;
      argv++;
      padding = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown option: " << argv[1] << endl;
      usage();
    }
  }

  // Create default interpolator
  if (interpolator == NULL) {
    interpolator = new irtkNearestNeighborInterpolateImageFunction;
  }

  // Isotropic?
  if(isotropic > 0) {
    // Resample image to isotropic voxels (smallest voxel dimension)
    double size;
    image->GetPixelSize(&xsize, &ysize, &zsize);
    size = xsize;
    size = (size < ysize) ? size : ysize;
    size = (size < zsize) ? size : zsize;
    xsize = size*isotropic; 
    ysize = size*isotropic; 
    zsize = size*isotropic;
    cerr << "Resampling image to isotropic voxel size (in mm): " << xsize << endl;
  } else {
    cout << "Resampling ... "; cout.flush();
  }
  switch (image->GetScalarType()) {
    case IRTK_VOXEL_UNSIGNED_CHAR: {
        if (padding == false) {
          irtkResampling<unsigned char> resampling(xsize, ysize, zsize);
          resampling.SetInput ((irtkGenericImage<unsigned char>*)(image));
          resampling.SetOutput((irtkGenericImage<unsigned char>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        } else {
          irtkResamplingWithPadding<unsigned char> resampling(xsize, ysize, zsize, padding_value);
          resampling.SetInput ((irtkGenericImage<unsigned char>*)(image));
          resampling.SetOutput((irtkGenericImage<unsigned char>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        }
        break;
      }
    case IRTK_VOXEL_SHORT: {
        if (padding == false) {
          irtkResampling<short> resampling(xsize, ysize, zsize);
          resampling.SetInput ((irtkGenericImage<short>*)(image));
          resampling.SetOutput((irtkGenericImage<short>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        } else {
          irtkResamplingWithPadding<short> resampling(xsize, ysize, zsize, padding_value);
          resampling.SetInput ((irtkGenericImage<short>*)(image));
          resampling.SetOutput((irtkGenericImage<short>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        }
        break;
      }
    case IRTK_VOXEL_UNSIGNED_SHORT: {
        if (padding == false) {
          irtkResampling<unsigned short> resampling(xsize, ysize, zsize);
          resampling.SetInput ((irtkGenericImage<unsigned short>*)(image));
          resampling.SetOutput((irtkGenericImage<unsigned short>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        } else {
          irtkResamplingWithPadding<unsigned short> resampling(xsize, ysize, zsize, padding_value);
          resampling.SetInput ((irtkGenericImage<unsigned short>*)(image));
          resampling.SetOutput((irtkGenericImage<unsigned short>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        }
        break;
      }
    case IRTK_VOXEL_INT: {
        if (padding == false) {
          irtkResampling<int> resampling(xsize, ysize, zsize);
          resampling.SetInput ((irtkGenericImage<int>*)(image));
          resampling.SetOutput((irtkGenericImage<int>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        } else {
          irtkResamplingWithPadding<int> resampling(xsize, ysize, zsize, padding_value);
          resampling.SetInput ((irtkGenericImage<int>*)(image));
          resampling.SetOutput((irtkGenericImage<int>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        }
        break;
      }
    case IRTK_VOXEL_FLOAT: {
        if (padding == false) {
          irtkResampling<float> resampling(xsize, ysize, zsize);
          resampling.SetInput ((irtkGenericImage<float>*)(image));
          resampling.SetOutput((irtkGenericImage<float>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        } else {
          irtkResamplingWithPadding<float> resampling(xsize, ysize, zsize, padding_value);
          resampling.SetInput ((irtkGenericImage<float>*)(image));
          resampling.SetOutput((irtkGenericImage<float>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        }
        break;
      }
    case IRTK_VOXEL_DOUBLE: {
        if (padding == false) {
          irtkResampling<double> resampling(xsize, ysize, zsize);
          resampling.SetInput ((irtkGenericImage<double>*)(image));
          resampling.SetOutput((irtkGenericImage<double>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        } else {
          irtkResamplingWithPadding<double> resampling(xsize, ysize, zsize, padding_value);
          resampling.SetInput ((irtkGenericImage<double>*)(image));
          resampling.SetOutput((irtkGenericImage<double>*)(image));
          resampling.SetInterpolator(interpolator);
          resampling.Run();
        }
        break;
      }
    default:
      cerr << "transformation: Unknown scalar type" << endl;
      exit(1);
  }
  ok = true;
  cerr << "done.."<<endl;

  // Save result
  image->Write(output_name);

  delete image;

  return 0;
}
