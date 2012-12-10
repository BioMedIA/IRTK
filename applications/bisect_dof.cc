/*=========================================================================
  bisect_dof.cc

  Library   : Image Registration Toolkit (IRTK)
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 23 Feb 2012 $
  Changes   : $Author: cl6311 $

=========================================================================*/

#define POINT_SCALING 100

#include <irtkRegistration.h>

char *dof_name_in = NULL;
char *dof_name_out = NULL;

void usage()
{
  cerr << "Usage: bisect_dof [dof_in] [dof_out] [-rigid/-affine]" << endl;
  cerr << "Calculates a square root of the transformation matrix." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, affine;

  irtkMatrix matrix;

  // Check arguments
  if (argc < 3) usage();

  // Parse arguments
  dof_name_in = argv[1];
  argc--;
  argv++;

  dof_name_out = argv[1];
  argc--;
  argv++;

  // Default
  affine = false;
  while (argc > 1) {
    ok = false;
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

  if (affine != true) {
	  // Create transformation;
	  irtkRigidTransformation transformation;
	  transformation.irtkTransformation::Read(dof_name_in);
	  transformation.Print();
	  irtkMatrix m = transformation.GetMatrix();
	  irtkMatrix msq = sqrtm(m);
	  transformation.PutMatrix(msq);

	  // Write transformation
	  transformation.irtkTransformation::Write(dof_name_out);
	  transformation.Print();
  }else{
	  // Create transformation;
	  irtkAffineTransformation transformation;
	  transformation.irtkTransformation::Read(dof_name_in);
	  transformation.Print();
	  irtkMatrix m = transformation.GetMatrix();
	  irtkMatrix msq = sqrtm(m);
	  transformation.PutMatrix(msq);

	  // Write transformation
	  transformation.irtkTransformation::Write(dof_name_out);
	  transformation.Print();
  }
  return 0;
}
