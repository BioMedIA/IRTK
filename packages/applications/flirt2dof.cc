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

#include <irtkFileToImage.h>

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

// Returns D*inv(M)*inv(D) in MJ's notes
irtkMatrix get_matrix(irtkGreyImage &image, int bReflect)
{
  int i, dims[3];
  double vdims[3];
  irtkMatrix C, D(4,4), D_inv(4,4), M(4,4);

  //
  image.GetPixelSize(&vdims[0], &vdims[1], &vdims[2]);
  dims[0] = image.GetX();
  dims[1] = image.GetY();
  dims[2] = image.GetZ();

  // D, and inverse
  for (i = 0; i < 3; i++) {
    D(i, i) = vdims[i];
  }
  D(3, 3) = 1;
  D_inv = D;
  D_inv.Invert();

  // M, and inverse
  for (i = 0; i < 4; i++) {
    M(i, i) = 1;
  }
  for (i = 0; i < 3; i++) {
    M(i, 3) = 0.5*(dims[i]-1);
  }
  M.Invert();

  // C
  C = D * M * D_inv;

  if (bReflect == true) {
    cerr << "Reflecting...\n";
    irtkMatrix Fy(4, 4);
    Fy(0, 0) = 1;
    Fy(1, 1) = -1;
    Fy(2, 2) = 1;
    Fy(3, 3) = 1;
    Fy(1, 3) = (image.GetY()-1)*vdims[1];
    C *= Fy;
  }

  cerr << "C: " << endl;
  C.Print();
  return C;
}

int main(int argc, char **argv)
{
  int i, j, ok, bTargetY, bSourceY;
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
  irtkMatrix matrix(4, 4), matrix_save(4, 4);
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      from >> val;
      matrix.Put(i, j, val);
    }
  }
  matrix_save = matrix;

  // Read in image pairs to deal with different coordinate systems
  irtkFileToImage *target_reader = irtkFileToImage::New(argv[1]);
  irtkGreyImage target(argv[1]);
  argc--;
  argv++;
  irtkFileToImage *source_reader = irtkFileToImage::New(argv[1]);
  irtkGreyImage source(argv[1]);
  argc--;
  argv++;

  // Set up affine transformation
  dofout_name = argv[1];
  argc--;
  argv++;

  // Check for Analyze images
  if (strcmp(target_reader->NameOfClass(), "irtkFileANALYZEToImage") == 0) {
    bTargetY = true;
  } else {
    bTargetY = false;
  }
  if (strcmp(source_reader->NameOfClass(), "irtkFileANALYZEToImage") == 0) {
    bSourceY = true;
  } else {
    bSourceY = false;
  }

  // Get conversion matrices
  irtkMatrix Ctgt = get_matrix(target, bTargetY);
  irtkMatrix Csrc = get_matrix(source, bSourceY);
  Ctgt.Invert();
  matrix.Invert();
  matrix = Csrc * matrix * Ctgt;

  // Set matrix
  irtkAffineTransformation transformation;
  transformation.PutMatrix(matrix); // updates parameters

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
      matrix_save.Print();
      matrix.Print();
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
