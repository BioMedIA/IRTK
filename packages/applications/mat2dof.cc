/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#define POINT_SCALING 100

#include <irtkRegistration.h>

char *mat_name = NULL;
char *dof_name = NULL;

void usage()
{
  cerr << "Usage: mat2dof [matfile] [doffile] [-rigid/-affine] [-invert] [-import]" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  double x, y, z, error;
  int i, ok, affine, import, invert;

  irtkMatrix matrix;
  irtkTransformation *transformation = NULL;
  irtkPointRegistration *registration = NULL;

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
  invert = False;
  affine = False;
  import = False;
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-import") == 0)) {
      argc--;
      argv++;
      import = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-rigid") == 0)) {
      argc--;
      argv++;
      affine = False;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-affine") == 0)) {
      argc--;
      argv++;
      affine = True;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read matrix
  if (import == True) {
    matrix.Import(mat_name, 4, 4);
    matrix.Print();
  } else {
    matrix.Read(mat_name);
    matrix.Print();
  }

  // Invert matrix if necessary
  if (invert == True) {
    matrix.Invert();
  }

  // Create homogeneous transformation from matrix
  irtkHomogeneousTransformation mattrans;
  mattrans.PutMatrix(matrix);

  // Create pointsets for registration
  irtkPointSet pset1;
  irtkPointSet pset2;
  for (x = -1; x <= 1; x += 2) {
    for (y = -1; y <= 1; y += 2) {
      for (z = -1; z <= 1; z += 2) {
        pset1.Add(irtkPoint(POINT_SCALING * x, POINT_SCALING * y, POINT_SCALING * z));
        pset2.Add(irtkPoint(POINT_SCALING * x, POINT_SCALING * y, POINT_SCALING * z));
      }
    }
  }

  // Transform point set
  mattrans.irtkTransformation::Transform(pset1);

  if (affine != True) {
    // Create transformation;
    transformation = new irtkRigidTransformation;
    // Create registration
    registration = new irtkPointRigidRegistration;
    // Setup registration
    irtkPointSet tmp1 = pset1;
    irtkPointSet tmp2 = pset2;
    registration->SetInput(&tmp1, &tmp2);
    registration->SetOutput(transformation);
    registration->Run();
  } else {
    // Create transformation
    transformation = new irtkAffineTransformation;
    // Create registration
    registration = new irtkPointAffineRegistration;
    // Setup registration
    irtkPointSet tmp1 = pset1;
    irtkPointSet tmp2 = pset2;
    registration->SetInput(&tmp1, &tmp2);
    registration->SetOutput(transformation);
    registration->Run();
  }

  // Calculate residual error
  transformation->irtkTransformation::Transform(pset2);

  error = 0;
  for (i = 0; i < pset1.Size(); i++) {
    irtkPoint p1 = pset1(i);
    irtkPoint p2 = pset2(i);
    error += sqrt(pow(double(p1._x - p2._x), 2.0) +
                  pow(double(p1._y - p2._y), 2.0) +
                  pow(double(p1._z - p2._z), 2.0));
  }
  error /= pset1.Size();
  if (error > 0.1) {
    cout << "RMS is " << error/pset1.Size() << " mm" << endl;
  }

  // Write transformation
  transformation->Write(dof_name);
  transformation->Print();
}
