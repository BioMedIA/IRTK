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

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataWriter.h>

#include <irtkTransformation.h>

#include <irtkSurfaceRegistration.h>

//some global variables
char *_target_name = NULL, *_source_name = NULL;
char *dofin_name = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "Usage: shreg [target] [source] [subdivisions] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-locator>          Locator: 0 = cell locator, 1 = point locator, 2 = kd-tree locator (default = 1)" << endl;
  cerr << "<-dofin name>       Name of input file" << endl;
  cerr << "<-dofout name>      Name of output file" << endl;
  cerr << "<-epsilon>          Value for espilon (default=0.01)" << endl;
  cerr << "<-symmetric>        Use symmetric distance (default OFF)" << endl;
  cerr << "<-ignoreedges>      Ignores edges in ICP (default OFF)" << endl;
  cerr << "<-iterations>       Number of 3D registration iterations (default 100)" << endl;
  cerr << "<-ds spacing>       Control point spacing" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, numberOfLevels, locatorType, iterations, ok;
  double epsilon, dx, dy, dz, start_box[3], end_box[3], *bounding_box;
  bool ignoreEdges, subdivide, symmetricDistance;

  if (argc < 3) {
    usage();
  }

  // Default parameters
  iterations = 100;
  locatorType = 1;
  epsilon = 0.01;
  ok = 0;
  ignoreEdges = False;
  subdivide = False;
  symmetricDistance = False;

  // Fix spacing
  dx = 20;
  dy = 20;
  dz = 20;

  // Parse filenames
  _target_name = argv[1];
  argv++;
  argc--;
  _source_name = argv[1];
  argv++;
  argc--;

  // Number of subdivisions
  numberOfLevels = atoi(argv[1]);
  argc--;
  argv++;

  // Parse remaining parameters
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-locator") == 0)) {
      argc--;
      argv++;
      locatorType = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-ignoreedges") == 0)) {
      argc--;
      argv++;
      ignoreEdges = true;
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
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-symmetric") == 0)) {
      argc--;
      argv++;
      symmetricDistance = True;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-epsilon") == 0)) {
      argc--;
      argv++;
      epsilon = atof(argv[1]);
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

  // Target pipeline
  cout << "Reading target ... " << _target_name;
  vtkPolyDataReader *target_reader = vtkPolyDataReader::New();
  target_reader->SetFileName(_target_name);
  target_reader->Modified();
  vtkPolyData *target = vtkPolyData::New();
  target = target_reader->GetOutput();
  target->Modified();
  target->Update();
  cout << "  (" << target->GetNumberOfPoints() << " points)" << endl;

  // Source pipeline
  cout << "Reading source ... " << _source_name;
  vtkPolyDataReader *source_reader = vtkPolyDataReader::New();
  source_reader->SetFileName(_source_name);
  source_reader->Modified();
  source_reader->Update();
  vtkPolyData *source = vtkPolyData::New();
  source = source_reader->GetOutput();
  source->Modified();
  cout << "  (" << source->GetNumberOfPoints() << " points)" << endl;

  if (ignoreEdges == True) {
    MarkBoundary(target);
    MarkBoundary(source);
  }

  // Create registration
  irtkSurfaceFreeFormRegistration *registration = new irtkSurfaceFreeFormRegistration;
  registration->SetInput(target, source);
  registration->SetNumberOfIterations(iterations);
  registration->SetEpsilon(epsilon);

  // Check if to do symmetric registration
  if (symmetricDistance == True) {
    // Create target locator
    irtkLocator *target_locator = new irtkLocator;
    target_locator->SelectLocatorType(locatorType);
    registration->SetTargetLocator(target_locator);
    // Create source locator
    irtkLocator *source_locator = new irtkLocator;
    source_locator->SelectLocatorType(locatorType);
    registration->SetSourceLocator(source_locator);
    registration->UseSymmetricDistance();
  } else {
    // Create source locator
    irtkLocator *source_locator = new irtkLocator;
    source_locator->SelectLocatorType(locatorType);
    registration->SetSourceLocator(source_locator);
  }

  // Create initial multi-level free-form deformation
  irtkMultiLevelFreeFormTransformation *mffd = NULL;

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
    // Otherwise use identity transformation to start
    mffd = new irtkMultiLevelFreeFormTransformation;
  }

  // Get bounding box of data
  bounding_box = target->GetBounds();
  start_box[0] = *(bounding_box)   - dx;
  start_box[1] = *(bounding_box+2) - dy;
  start_box[2] = *(bounding_box+4) - dz;
  end_box[0]   = *(bounding_box+1) + dx;
  end_box[1]   = *(bounding_box+3) + dy;
  end_box[2]   = *(bounding_box+5) + dz;

  // Create transformation
  double xaxis[3] = {1, 0, 0};
  double yaxis[3] = {0, 1, 0};
  double zaxis[3] = {0, 0, 1};
  irtkBSplineFreeFormTransformation *affd =  new
  irtkBSplineFreeFormTransformation(start_box[0], start_box[1], start_box[2],
                                    end_box[0], end_box[1], end_box[2], dx, dy, dz, xaxis, yaxis, zaxis);

  // Add ffd
  mffd->PushLocalTransformation(affd);

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
      irtkBSplineFreeFormTransformation(start_box[0], start_box[1], start_box[2],
                                        end_box[0], end_box[1], end_box[2], dx, dy, dz, xaxis, yaxis, zaxis);
      mffd->PushLocalTransformation(affd);
    } else {
      // Extract current transformation level
      irtkBSplineFreeFormTransformation *affd =
        (irtkBSplineFreeFormTransformation *)mffd->GetLocalTransformation(0);

      // Subdivide
      affd->Subdivide();
    }
  }

  // Set up registration and run
  registration->SetOutput(mffd);
  registration->Run();

  // Write rigid transformation
  if (dofout_name != NULL) {
    mffd->irtkTransformation::Write(dofout_name);
  }
}


#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
