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

#include <irtkTransformation.h>

// Default filenames
char *target_name = NULL, *source_name = NULL, *output_name = NULL;

void usage()
{
  cerr << "Usage: substract [target] [source] [output] <options>" << endl;
  cerr << "where <options> is one or more of the following:\n"    << endl;
  cerr << "<-Tp value>          Padding value in target" << endl;
  cerr << "<-Sp value>          Padding value in source" << endl;
  cerr << "<-no_norm>           Don't normalize images" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, ok, n, no_norm, target_padding, source_padding;
  double a, b, cov, var, x_avg, y_avg;

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse image
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  // Fix padding
  target_padding = MIN_GREY;
  source_padding = MIN_GREY;
  no_norm = False;
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      target_padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sp") == 0)) {
      argc--;
      argv++;
      source_padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-no_norm") == 0)) {
      argc--;
      argv++;
      no_norm = True;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can't parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read image
  cout << "Reading image ... "; cout.flush();
  irtkGreyImage target(target_name);
  cout << "done" << endl;

  cout << "Reading image ... "; cout.flush();
  irtkGreyImage source(source_name);
  cout << "done" << endl;

  // Check whether images are of same size
  if (!((target.GetX() == source.GetX()) &&
        (target.GetY() == source.GetY()) &&
        (target.GetZ() == source.GetZ()))) {
    cerr << "subtract: Images must have same size" << endl;
    exit(1);
  }

  if (no_norm == False) {

    // Normalize and subtract
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
          if ((source(i, j, k) > source_padding) && (target(i, j, k) > target_padding)) {
            target(i, j, k) = target(i, j, k) - round(a + b * source(i, j, k));
          } else {
            target(i, j, k) = 0;
          }
        }
      }
    }
  } else {

    // Subtract images
    for (k = 0; k < source.GetZ(); k++) {
      for (j = 0; j < source.GetY(); j++) {
        for (i = 0; i < source.GetX(); i++) {
          if ((source(i, j, k) > source_padding) && (target(i, j, k) > target_padding)) {
            target(i, j, k) = target(i, j, k) - source(i, j, k);
          } else {
            target(i, j, k) = 0;
          }
        }
      }
    }
  }

  // Write the substracted images
  target.Write(output_name);

  // Print out some info
  short min, max;
  target.GetMinMax(&min, &max);
  target.Print();
  cerr << "MinMax: " << min << " " << max << endl;
}

