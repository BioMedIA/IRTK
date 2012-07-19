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
#include <irtkHistogram.h>
#include <irtkSegmentationFunction.h>
#include <irtkGaussianBlurring.h>
#include <irtkNormalizeNyul.h>

// Default filenames
char *target_name = NULL, *source_name = NULL, *sourceoutput_name = NULL, *targetoutput_name = NULL;

void usage()
{
  cerr << "Usage: normalize [reference image] [image to be normalized] [output] ";
  cerr << "<options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-Tp value>        Target padding value" << endl;
  cerr << "<-Sp value>        Source padding value" << endl;
  cerr << "<-equalize paddingvalue targetoutput>  equalize the image's histogram and backproject" << endl;
  cerr << "<-piecewise>       use a piecewise linear function as suggested by Nyul et al." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, n, ok;
  int target_padding, source_padding;
  double a, b, cov, var, x_avg, y_avg,x,y,z;
  bool piecewise = false;

  if(argc < 4) usage();

  // Parse image
  target_name  = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;
  sourceoutput_name = argv[1];
  argc--;
  argv++;

  source_padding = MIN_GREY;
  target_padding = MIN_GREY;

  // Read image
  cout << "Reading target image ... "; cout.flush();
  irtkRealImage target (target_name);
  cout << "done" << endl;

  // Read image
  cout << "Reading source image ... "; cout.flush();
  irtkRealImage source(source_name);
  cout << "done" << endl;

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      target_padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sp") == 0)) {
      argc--;
      argv++;
      source_padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-equalize") == 0)) {
		argc--;
		argv++;
		double padding = atoi(argv[1]);
        if(padding < 0){source = nn.GetOutput();
            cout << "equalize padding value < 0 are you sure about it?" << endl;
        }
		argc--;
		argv++;
		targetoutput_name = argv[1];
		argc--;
		argv++;
		double min,max;
		target.GetMinMaxAsDouble(&min,&max);
        if (min < padding) min = padding;
		irtkImageHistogram_1D<irtkRealPixel> histogram;
		histogram.Evaluate(&target,padding);
		histogram.Equalize(min,max);
		histogram.BackProject(&target);
        source.GetMinMaxAsDouble(&min,&max);
        if (min < padding) min = padding;
		histogram.Evaluate(&source,padding);
		histogram.Equalize(min,max);
		histogram.BackProject(&source);

		ok = true;
	}
    if ((ok == false) && (strcmp(argv[1], "-piecewise") == 0)) {
      argc--;
      argv++;
      piecewise = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }
  if(piecewise == false){
	  n = 0;
	  x_avg = 0;
	  y_avg = 0;
	  for (k = 0; k < source.GetZ(); k++) {
		for (j = 0; j < source.GetY(); j++) {
		  for (i = 0; i < source.GetX(); i++) {
			  x = i; y = j; z = k;
			  source.ImageToWorld(x,y,z);
			  target.WorldToImage(x,y,z);
			  x = round(x); y = round(y); z = round(z);
			  if(x >= 0 && x < target.GetX()
				  && y >= 0 && y < target.GetY()
				  && z >= 0 && z < target.GetZ()){
					  if ((source(i, j, k) > source_padding) && (target(x,y,z) > target_padding)) {
						  n++;
						  x_avg += source(i, j, k);
						  y_avg += target(x, y, z);
					  }
			  }
		  }
		}
	  }
	  if (n == 0) {
		cerr << "Normalize: Number of samples should be larger than zero" << endl;
		exit(1);
	  }

	  cov = 0;
	  var = 0;
	  for (k = 0; k < source.GetZ(); k++) {
		  for (j = 0; j < source.GetY(); j++) {
			  for (i = 0; i < source.GetX(); i++) {
				  x = i; y = j; z = k;
				  source.ImageToWorld(x,y,z);
				  target.WorldToImage(x,y,z);
				  x = round(x); y = round(y); z = round(z);
				  if(x >= 0 && x < target.GetX()
					  && y >= 0 && y < target.GetY()
					  && z >= 0 && z < target.GetZ()){
						  if ((source(i, j, k) > source_padding) && (target(x, y, z) > target_padding)) {
							  cov += (source(i, j, k) - x_avg) * (target(x, y, z) - y_avg);
							  var += (source(i, j, k) - x_avg) * (source(i, j, k) - x_avg);
						  }
				  }
			  }
		  }
	  }
	  cov /= n;
	  var /= n;
	  b = cov / var;
	  a = y_avg - b * x_avg;

	  cout << "Scaling = " << b << endl;
	  cout << "Offset  = " << a << endl;

	  for (k = 0; k < source.GetZ(); k++) {
		for (j = 0; j < source.GetY(); j++) {
		  for (i = 0; i < source.GetX(); i++) {
			if (source(i, j, k) > source_padding) {
				source(i, j, k) = round(a + b * source(i, j, k));
			}
		  }
		}
	  }
  }
  else{ //if piecewise == true -> use piecewise linear approach proposed by Nyul et al.
	  irtkNormalizeNyul nn(source, target);
	  nn.SetPadding(source_padding, target_padding);
	  nn.Run();
	  source = nn.GetOutput();
  }

  // Write image
  source.Write(sourceoutput_name);
  if(targetoutput_name){
	target.Write(targetoutput_name);
  }
}
