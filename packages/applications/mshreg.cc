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

#include <irtkMultipleSurfaceRegistration.h>

//some global variables
char *_target_name = NULL, *_source_name = NULL;
char *dofin_name = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "Usage: mshreg [no. of surfaces] [target] [source] [subdivisions] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-locator>          Locator: 0 = cell locator, 1 = point locator, 2 = kd-tree locator (default = 1)" << endl;
  cerr << "<-dofin name>       Name of input file" << endl;
  cerr << "<-dofout name>      Name of output file" << endl;
  cerr << "<-epsilon>          Value for espilon (default=0.01)" << endl;
  cerr << "<-symmetric>        Use symmetric distance (default OFF)" << endl;
  cerr << "<-iterations>       Number of 3D registration iterations (default 100)" << endl;
  cerr << "<-ds spacing>       Control point spacing" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, numberOfLevels, no_surfaces, locatorType, iterations, ok;
  double epsilon, dx, dy, dz, start_box[3], end_box[3], *bounding_box;
  bool subdivide, symmetricDistance;

  if (argc < 3) {
    usage();
  }

  // Default parameters
  iterations = 100;
  locatorType = 1;
  epsilon = 0.01;
  ok = 0;
  subdivide = false;
  symmetricDistance = false;

  // Fix spacing
  dx = 20;
  dy = 20;
  dz = 20;

  // Parse number of surfaces
  no_surfaces = atoi(argv[1]);
  argv++;
  argc--;

  vtkPolyDataReader *target_reader = vtkPolyDataReader::New();
  vtkCleanPolyData *target_cleaner = vtkCleanPolyData::New();
  vtkPolyData **target = new vtkPolyData*[no_surfaces];
  for (i = 0; i < no_surfaces; i++) {
    // Parse target filename
    _target_name = argv[1];
    argv++;
    argc--;
    // Target pipeline
    cout << "Reading target ... " << _target_name << endl;
    target_reader->SetFileName(_target_name);
    target_reader->Modified();
    target_cleaner->SetInput(target_reader->GetOutput());
    target_cleaner->Modified();
    target_cleaner->Update();
    target[i] = target_cleaner->GetOutput();
    target[i]->Update();
  }

  vtkPolyDataReader *source_reader = vtkPolyDataReader::New();
  vtkCleanPolyData *source_cleaner = vtkCleanPolyData::New();
  vtkPolyData **source = new vtkPolyData*[no_surfaces];
  for (i = 0; i < no_surfaces; i++) {
    // Parse source filename
    _source_name = argv[1];
    argv++;
    argc--;
    // Source pipeline
    cout << "Reading source ... " << _source_name << endl;
    source_reader->SetFileName(_source_name);
    source_reader->Modified();
    source_cleaner->SetInput(source_reader->GetOutput());
    source_cleaner->Modified();
    source_cleaner->Update();
    source[i] = source_cleaner->GetOutput();
    source[i]->Update();
  }

  // Number of subdivisions
  numberOfLevels = atoi(argv[1]);
  argc--;
  argv++;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-locator") == 0)) {
      argc--;
      argv++;
      locatorType = atoi(argv[1]);
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
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-symmetric") == 0)) {
      argc--;
      argv++;
      symmetricDistance = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-epsilon") == 0)) {
      argc--;
      argv++;
      epsilon = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-subdivide") == 0)) {
      argc--;
      argv++;
      subdivide = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ds") == 0)) {
      argc--;
      argv++;
      dx = atof(argv[1]);
      dy = atof(argv[1]);
      dz = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkMultipleSurfaceFreeFormRegistration *registration = new irtkMultipleSurfaceFreeFormRegistration;
  registration->SetInput(target, source);
  registration->SetNumberOfIterations(iterations);
  registration->SetNumberOfSurfaces(no_surfaces);
  registration->SetEpsilon(epsilon);

  // Check if to do symmetric registration
  if (symmetricDistance == true) {
    // Create target and source locators
    irtkLocator **target_locator = new irtkLocator *[no_surfaces];
    irtkLocator **source_locator = new irtkLocator *[no_surfaces];
    for (i = 0; i < no_surfaces; i++) {
      target_locator[i] = new irtkLocator;
      target_locator[i]->SelectLocatorType(locatorType);
      source_locator[i] = new irtkLocator;
      source_locator[i]->SelectLocatorType(locatorType);
    }
    registration->SetTargetLocator(target_locator);
    registration->SetSourceLocator(source_locator);
    registration->UseSymmetricDistance();
  } else {
    // Create source locators
    irtkLocator **source_locator = new irtkLocator *[no_surfaces];
    for (i = 0; i < no_surfaces; i++) {
      source_locator[i] = new irtkLocator;
      source_locator[i]->SelectLocatorType(locatorType);
    }
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
  bounding_box = target[0]->GetBounds();
  start_box[0] = *(bounding_box)   - dx;
  start_box[1] = *(bounding_box+2) - dy;
  start_box[2] = *(bounding_box+4) - dz;
  end_box[0]   = *(bounding_box+1) + dx;
  end_box[1]   = *(bounding_box+3) + dy;
  end_box[2]   = *(bounding_box+5) + dz;
  for (i = 1; i < no_surfaces; i++) {
    bounding_box = target[i]->GetBounds();
    if ((*(bounding_box)   - dx) < start_box[0]) start_box[0] = *(bounding_box)   - dx;
    if ((*(bounding_box+2) - dy) < start_box[1]) start_box[1] = *(bounding_box+2) - dy;
    if ((*(bounding_box+4) - dz) < start_box[2]) start_box[2] = *(bounding_box+4) - dy;
    if ((*(bounding_box+1) + dx) > end_box[0])   end_box[0]   = *(bounding_box+1) + dx;
    if ((*(bounding_box+3) + dy) > end_box[1])   end_box[1]   = *(bounding_box+3) + dy;
    if ((*(bounding_box+5) + dz) > end_box[2])   end_box[2]   = *(bounding_box+5) + dz;
  }

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

    if (subdivide == false) {
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
