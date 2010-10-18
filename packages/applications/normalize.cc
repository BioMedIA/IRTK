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

// Default filenames
char *target_name = NULL, *source_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: normalize [reference image] [image to be normalized] [output] ";
  cerr << "<options>\n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-Tp value>        Target padding value" << endl;
  cerr << "<-Sp value>        Source padding value" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, n, ok;
  int target_padding, source_padding;
  double a, b, cov, var, x_avg, y_avg;

  // Parse image
  target_name  = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  source_padding = 0;
  target_padding = MIN_GREY;

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
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read image
  cout << "Reading target image ... "; cout.flush();
  irtkGreyImage target (target_name);
  cout << "done" << endl;

  // Read image
  cout << "Reading source image ... "; cout.flush();
  irtkGreyImage source(source_name);
  cout << "done" << endl;

  n = 0;
  x_avg = 0;
  y_avg = 0;
  for (k = 0; k < source.GetZ(); k++) {
    for (j = 0; j < source.GetY(); j++) {
      for (i = 0; i < source.GetX(); i++) {
        if ((source(i, j, k) > source_padding) && (target(i, j, k) > target_padding)) {
          n++;
          x_avg += source(i, j, k);
          y_avg += target(i, j, k);
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
        if ((source(i, j, k) > source_padding) && (target(i, j, k) > target_padding)) {
          cov += (source(i, j, k) - x_avg) * (target(i, j, k) - y_avg);
          var += (source(i, j, k) - x_avg) * (source(i, j, k) - x_avg);
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

  // Write image
  source.Write(output_name);
}
