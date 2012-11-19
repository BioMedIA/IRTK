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

#include <irtkTransformation.h>

char *dofout_name = NULL;

void usage()
{
  cerr << "Usage: flirt2dof [matrix_in.raw] [target] [source] [dofout] <options>\n";
  cerr << "where <options> can be one or more of the following:\n";
  cerr << "<-info>          \t\t Print transformation info\n";
  cerr << "<-print>         \t\t Print matrix to screen\n\n";
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, ok;
  double val;

  // Check command line
  if (argc < 5) {
    usage();
  }

  // Read in matrix
  ifstream from(argv[1]);
  if (!from) {
    cerr << "Can't open file " << argv[1] << "\n";
    exit(1);
  }
  argc--;
  argv++;

  irtkMatrix flirt_matrix(4, 4);
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      from >> val;
      flirt_matrix.Put(i, j, val);
    }
  }


  // Read in image pairs to deal with different coordinate systems
  irtkGreyImage target(argv[1]);
  argc--;
  argv++;

  irtkGreyImage source(argv[1]);
  argc--;
  argv++;

  // Set up affine transformation
  dofout_name = argv[1];
  argc--;
  argv++;

  irtkMatrix w2iTgt = target.GetWorldToImageMatrix();

  irtkMatrix i2wSrc = source.GetImageToWorldMatrix();

  irtkMatrix samplingTarget(4,4);
  irtkMatrix samplingSource(4,4);

  samplingTarget.Ident();
  samplingSource.Ident();

  samplingTarget(0, 0) = target.GetXSize();
  samplingTarget(1, 1) = target.GetYSize();
  samplingTarget(2, 2) = target.GetZSize();

  samplingSource(0, 0) = source.GetXSize();
  samplingSource(1, 1) = source.GetYSize();
  samplingSource(2, 2) = source.GetZSize();

  irtkMatrix outputMatrix(4, 4);

  // See fsl/src/newimage/newimagefns.cc : raw_affine_transform(..)
  flirt_matrix.Invert();
  samplingSource.Invert();
  outputMatrix = i2wSrc * samplingSource * flirt_matrix * samplingTarget * w2iTgt;

  // Set matrix
  irtkAffineTransformation transformation;
  transformation.PutMatrix(outputMatrix); // updates parameters

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-info") == 0)) {
      transformation.Print();
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-print") == 0)) {
      // Print flirt matrix that was read in.
      flirt_matrix.Invert();
      flirt_matrix.Print();
      outputMatrix.Print();
      argc--;
      argv++;
      ok = true;
    }
    if (!ok) {
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Write out
  transformation.irtkTransformation::Write(dofout_name);

  return 0;
}
