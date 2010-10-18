/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifdef HAS_VTK

#include <irtkRegistration.h>

// Default filenames
char *source_name = NULL, *target_name = NULL;
char *dofin_name  = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "Usage: pareg [target] [source]\n" << endl;
  cerr << "<-dofout    file>    Final transformation estimate" << endl;
  cerr << "<-dofin     file>    Start transformation estimate" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, i;
  double error;
  irtkPointSet target, source;

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Parse source and target point lists
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  // Read point set
  target.ReadVTK(target_name);

  // Read point set
  source.ReadVTK(source_name);

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofout") == 0)) {
      argc--;
      argv++;
      dofout_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Create registration filter
  irtkPointAffineRegistration *registration = new irtkPointAffineRegistration;

  // Create initial affine transformation
  irtkAffineTransformation *transformation = NULL;

  if (dofin_name != NULL) {
    irtkTransformation *transform = irtkTransformation::New(dofin_name);
    if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
      transformation = new irtkAffineTransformation(*((irtkRigidTransformation *)transform));
    } else {
      if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
        transformation = new irtkAffineTransformation(*((irtkAffineTransformation *)transform));
      } else {
        cerr << "Input transformation is not of type rigid or affine " << endl;
        exit(1);
      }
    }
    delete transform;
  } else {
    // Otherwise use identity transformation to start
    transformation = new irtkAffineTransformation;
  }

  // Set input and output for the registration filter
  irtkPointSet tmp1 = target;
  irtkPointSet tmp2 = source;
  registration->SetInput(&tmp1, &tmp2);
  registration->SetOutput(transformation);

  // Run registration filter
  registration->Run();

  // Calculate residual error
  transformation->irtkTransformation::Transform(target);

  error = 0;
  for (i = 0; i < target.Size(); i++) {
    irtkPoint p1 = target(i);
    irtkPoint p2 = source(i);
    error += sqrt(pow(double(p1._x - p2._x), 2.0) +
                  pow(double(p1._y - p2._y), 2.0) +
                  pow(double(p1._z - p2._z), 2.0));
  }
  cout << "RMS is " << error/target.Size() << " mm" << endl;

  // Write the final transformation estimate
  if (dofout_name != NULL) transformation->irtkTransformation::Write(dofout_name);
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

