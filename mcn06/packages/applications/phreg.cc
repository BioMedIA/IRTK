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
  cerr << "Usage: phreg [target] [source] [subdivisions] <options>\n" << endl;
  cerr << "<-dofout file>       Final transformation estimate" << endl;
  cerr << "<-dofin  file>       Start transformation estimate" << endl;
  cerr << "<-ds spacing>        Control point spacing" << endl;
  cerr << "<-subdivide>         Subdivide" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, i, numberOfLevels, subdivide;
  double error, dx, dy, dz;
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

  // Number of subdivisions
  numberOfLevels = atoi(argv[1]);
  argc--;
  argv++;

  // Default options
  subdivide = False;

  // Fix spacing
  dx = 20;
  dy = 20;
  dz = 20;

  // Read point set
  target.ReadVTK(target_name);

  // Read point set
  source.ReadVTK(source_name);

  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dofout") == 0)) {
      argc--;
      argv++;
      dofout_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-subdivide") == 0)) {
      argc--;
      argv++;
      subdivide = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-ds") == 0)) {
      argc--;
      argv++;
      dx = atof(argv[1]);
      dy = atof(argv[1]);
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Create initial multi-level free-form deformation
  irtkMultiLevelFreeFormTransformation *mffd = NULL;

  // Read transformation
  if (dofin_name != NULL) {
    irtkTransformation *transform = irtkTransformation::New(dofin_name);
    if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
      mffd = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
    } else {
      if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
        mffd = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
      } else {
        if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
          mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
        } else {
          cerr << "Input transformation is not of type rigid, affine " << endl;
          cerr << "or multi-level free form deformation" << endl;
          exit(1);
        }
      }
    }
    delete transform;
  } else {
    mffd = new irtkMultiLevelFreeFormTransformation;
  }

  // Get bounding box of data
  irtkPoint p1, p2;
  target.BoundingBox(p1, p2);
  p1 -= irtkPoint(dx, dy, dz);
  p2 += irtkPoint(dx, dy, dz);

  // Create transformation
  double xaxis[3] = {1, 0, 0};
  double yaxis[3] = {0, 1, 0};
  double zaxis[3] = {0, 0, 1};
  irtkBSplineFreeFormTransformation *affd =  new
  irtkBSplineFreeFormTransformation(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z, dx, dy, dz, xaxis, yaxis, zaxis);

  // Add ffd
  mffd->PushLocalTransformation(affd);

  // Create registration
  irtkPointFreeFormRegistration *registration = new irtkPointFreeFormRegistration;
  registration->SetInput(&target, &source);

  for (i = 0; i < numberOfLevels-1; i++) {

    // Set up registration and run
    registration->SetOutput(mffd);
    registration->Run();

    if (subdivide == False) {
      // Add transformation
      dx = dx/2.0;
      dy = dy/2.0;
      dz = dz/2.0;
      irtkBSplineFreeFormTransformation *affd =  new
      irtkBSplineFreeFormTransformation(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z, dx, dy, dz, xaxis, yaxis, zaxis);

      mffd->PushLocalTransformation(affd);
    } else {
      // Extract current transformation level
      irtkFreeFormTransformation *affd = mffd->GetLocalTransformation(0);

      // Subdivide
      affd->Subdivide();
    }
  }

  // Run registration filter
  registration->SetOutput(mffd);
  registration->Run();

  // Calculate residual error
  mffd->irtkTransformation::Transform(target);

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
  if (dofout_name != NULL) mffd->irtkTransformation::Write(dofout_name);
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

