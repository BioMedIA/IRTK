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

char *mat_name = NULL;
char *dof_name = NULL;

void usage()
{
  cerr << "Usage: mat2dof [matfile] [doffile] [-rigid/-affine] [-invert] [-import]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, affine, import, invert;

  irtkMatrix matrix;

  // Check arguments
  if (argc < 3) usage();

  // Parse arguments
  mat_name  = argv[1];
  argc--;
  argv++;
  dof_name = argv[1];
  argc--;
  argv++;

  // Default
  invert = false;
  affine = true;
  import = false;
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-import") == 0)) {
      argc--;
      argv++;
      import = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-rigid") == 0)) {
      argc--;
      argv++;
      affine = false;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-affine") == 0)) {
      argc--;
      argv++;
      affine = true;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read matrix
  if (import == true) {
    matrix.Import(mat_name, 4, 4);
    matrix.Print();
  } else {
    matrix.Read(mat_name);
    matrix.Print();
  }

  // Invert matrix if necessary
  if (invert == true) {
    matrix.Invert();
  }

  if (affine != true) {
	  // Create transformation;
	  irtkRigidTransformation transformation;
	  transformation.PutMatrix(matrix);
	  transformation.UpdateParameter();
	  // Write transformation
	  transformation.irtkTransformation::Write(dof_name);
	  transformation.Print();
  }else{
	  irtkAffineTransformation transformation;
	  transformation.PutMatrix(matrix);
	  transformation.UpdateParameter();
	  // Write transformation
	  transformation.irtkTransformation::Write(dof_name);
	  transformation.Print();
  }

}
