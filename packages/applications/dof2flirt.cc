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

char *matrix_outname = NULL;

void usage()
{
  cerr << "Usage: dof2flirt [dofin] [target] [source] [matrix_out.raw] <options>\n";
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

  if (bReflect == True) {
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
  irtkMatrix matrix(4, 4);

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
  irtkFileToImage<irtkGreyPixel> *target_reader =
    irtkFileToImage<irtkGreyPixel>::New(argv[1]);
  irtkGreyImage target(argv[1]);
  argc--;
  argv++;
  irtkFileToImage<irtkGreyPixel> *source_reader =
    irtkFileToImage<irtkGreyPixel>::New(argv[1]);
  irtkGreyImage source(argv[1]);
  argc--;
  argv++;

  // Read in output name
  matrix_outname = argv[1];
  argc--;
  argv++;

  // Check for Analyze images
  if (strcmp(target_reader->NameOfClass(), "irtkFileANALYZEToImage") == 0) {
    bTargetY = True;
  } else {
    bTargetY = False;
  }
  if (strcmp(source_reader->NameOfClass(), "irtkFileANALYZEToImage") == 0) {
    bSourceY = True;
  } else {
    bSourceY = False;
  }

  // Get affine matrix
  matrix = transformation.GetMatrix();

  // Get conversion matrices
  irtkMatrix Ctgt = get_matrix(target, bTargetY);
  irtkMatrix Csrc = get_matrix(source, bSourceY);
  Ctgt.Invert();
  matrix = Ctgt * matrix * Csrc;

  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-info") == 0)) {
      transformation.Print();
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-print") == 0)) {
      matrix.Print();
      argc--;
      argv++;
      ok = True;
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
      to << matrix.Get(i, j) << " ";
    }
    to << endl;
  }

  return 0;
}
