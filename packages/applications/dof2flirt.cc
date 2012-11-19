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

char *matrix_outname = NULL;

void usage()
{
  cerr << "Usage: dof2flirt [dofin] [target] [source] [matrix_out.raw] <options>\n";
  cerr << "where <options> can be one or more of the following:\n";
  cerr << "<-info>          \t\t Print transformation info\n";
  cerr << "<-print>         \t\t Print matrix to screen\n\n";
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, ok;
  irtkMatrix dof_matrix(4, 4);

  // Check command line
  if (argc < 5) {
    usage();
  }

  // Read input transformation
  irtkAffineTransformation transformation;
  transformation.irtkTransformation::Read(argv[1]);
  argc--;
  argv++;

  // Read in image pairs to deal with different coordinate systems
  irtkGreyImage target(argv[1]);
  argc--;
  argv++;

  irtkGreyImage source(argv[1]);
  argc--;
  argv++;

  // Read in output name
  matrix_outname = argv[1];
  argc--;
  argv++;


  // Get affine matrix
  dof_matrix = transformation.GetMatrix();

  irtkMatrix i2wTgt = target.GetImageToWorldMatrix();

  irtkMatrix w2iSrc = source.GetWorldToImageMatrix();

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

  irtkMatrix flirt_matrix(4, 4);

  // See fsl/src/newimage/newimagefns.cc : raw_affine_transform(..)
  samplingTarget.Invert();
  flirt_matrix = samplingSource * w2iSrc * dof_matrix * i2wTgt * samplingTarget;
  flirt_matrix.Invert();


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
      flirt_matrix.Print();
      argc--;
      argv++;
      ok = true;
    }
    if (!ok) {
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Write out  matrix
  ofstream to(matrix_outname);
  if (!to) {
    cerr << "Can't open file " << matrix_outname << "\n";
    exit(1);
  }
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      to << flirt_matrix.Get(i, j) << " ";
    }
    to << endl;
  }

  return 0;
}
