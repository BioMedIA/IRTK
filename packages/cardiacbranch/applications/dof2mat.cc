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

char *matout_name = NULL;

void usage()
{
  cerr << "Usage: dof2mat [doffile] [-matout matrixfile] [-invert]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  irtkMatrix matrix;
  irtkAffineTransformation transformation;

  if (argc < 2) {
    usage();
  }

  // Read transformation
  transformation.irtkTransformation::Read(argv[1]);
  argc--;
  argv++;

  // Convert to matrix
  matrix = transformation.GetMatrix();

// Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-matout") == 0)) {
      argc--;
      argv++;
      matout_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      matrix.Invert();
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Print matrix
  matrix.Print();

  // Write matrix if necessary
  if (matout_name != NULL) {
    matrix.Write(matout_name);
  }
}
