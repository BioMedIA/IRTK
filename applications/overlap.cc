/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkImageFunction.h>

#include <irtkHistogram.h>

#include <irtkTransformation.h>

// Default filenames
char *source_name = NULL, *target_name = NULL;

void usage()
{
  cerr << "Usage: overlap [target] [source] <options>\n" << endl;
  cerr << "where <options> is one or more of the following: \n" << endl;
  cerr << "<-Rx1 pixel>       Region of interest" << endl;
  cerr << "<-Ry1 pixel>       Region of interest" << endl;
  cerr << "<-Rz1 pixel>       Region of interest" << endl;
  cerr << "<-Rx2 pixel>       Region of interest" << endl;
  cerr << "<-Ry2 pixel>       Region of interest" << endl;
  cerr << "<-Rz2 pixel>       Region of interest" << endl;
  cerr << "<-Tp  value>       Padding value in target" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, i1, i2, j1, j2, k1, k2, ok, max;
  irtkGreyPixel padding;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  // Read target and source image
  irtkGreyImage target(target_name);
  irtkGreyImage source(source_name);

  // Default padding
  padding = MIN_GREY;

  // Fix ROI
  i1 = 0;
  j1 = 0;
  k1 = 0;
  i2 = target.GetX();
  j2 = target.GetY();
  k2 = target.GetZ();

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      padding = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      i1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      i2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      j1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      j2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      k1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      k2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // If there is an region of interest, use it
  if ((i1 != 0) || (i2 != target.GetX()) ||
      (j1 != 0) || (j2 != target.GetY()) ||
      (k1 != 0) || (k2 != target.GetZ())) {
    target = target.GetRegion(i1, j1, k1, i2, j2, k2);
    source = source.GetRegion(i1, j1, k1, i2, j2, k2);
  }

  // Compute maximum value
  max = 0;
  for (k = 0; k < target.GetZ(); k++) {
    for (j = 0; j < target.GetY(); j++) {
      for (i = 0; i < target.GetX(); i++) {
        if (target(i, j, k) > max) max = target(i, j, k);
        if (source(i, j, k) > max) max = source(i, j, k);
      }
    }
  }

  irtkHistogram_1D<int> histogramA(max);
  irtkHistogram_1D<int> histogramB(max);
  irtkHistogram_1D<int> histogramAB(max);

  // Compute overlap
  for (k = 0; k < target.GetZ(); k++) {
    for (j = 0; j < target.GetY(); j++) {
      for (i = 0; i < target.GetX(); i++) {
        histogramA.Add(target(i, j, k));
        histogramB.Add(source(i, j, k));
        if (target(i, j, k) == source(i, j, k)) histogramAB.Add(target(i, j, k));
      }
    }
  }

  double w, si = 0;
  for (i = 1; i < histogramA.NumberOfBins(); i++) {
    w = histogramA(i) / double(histogramA.NumberOfSamples() - histogramA(0));
    si += w * 2.0 * histogramAB(i) / double(histogramA(i) + histogramB(i));
  }

  // Print average SI
  cout << "Average SI = " << si << endl;
}

