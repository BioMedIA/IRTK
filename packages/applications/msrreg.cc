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

char *_target_name = NULL, *_source_name = NULL;
char *dofin_name = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "Usage: msrreg [number of surfaces] [target] [source] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-locator>          Locator: 0 = cell locator, 1 = point locator, 2 = kd-tree locator (default = 1)" << endl;
  cerr << "<-optimizer>        Optimizer: 0 = gradient descent, 1 = conjugate gradient, 2 = downhill (default = 0)" << endl;
  cerr << "<-dofin name>       Name of input file" << endl;
  cerr << "<-dofout name>      Name of output file" << endl;
  cerr << "<-epsilon>          Value for espilon (default=0.01)" << endl;
  cerr << "<-clean>            Clean polydata (default OFF)" << endl;
  cerr << "<-symmetric>        Use symmetric distance (default OFF)" << endl;
  cerr << "<-invert>           Save the inverse transformation (default OFF)" << endl;
  cerr << "<-iterations>       Number of 3D registration iterations (default 100)" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  float epsilon;
  int i, locatorType, iterations, no_surfaces, ok;
  bool clean, invert, symmetricDistance;

  if (argc < 3) {
    usage();
  }

  // Default parameters
  iterations = 100;
  locatorType = 1;
  epsilon = 0.01;
  ok = 0;
  clean = false;
  invert = false;
  symmetricDistance = false;

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
    if (clean == true) {
      target_cleaner->SetInput(target_reader->GetOutput());
      target_cleaner->Modified();
      target_cleaner->Update();
      target[i] = target_cleaner->GetOutput();
    } else {
      target[i] = target_reader->GetOutput();
    }
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
    if (clean == true) {
      source_cleaner->SetInput(source_reader->GetOutput());
      source_cleaner->Modified();
      source_cleaner->Update();
      source[i] = source_cleaner->GetOutput();
    } else {
      source[i] = source_reader->GetOutput();
    }
    source[i]->Update();
  }

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
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
    if ((ok == false) && (strcmp(argv[1], "-clean") == 0)) {
      argc--;
      argv++;
      clean = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = true;
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
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--;
      argv++;
      iterations = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Create registration
  irtkMultipleSurfaceRigidRegistration *registration = new  irtkMultipleSurfaceRigidRegistration;
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

  // Create rigid transformation
  irtkRigidTransformation *transformation = new irtkRigidTransformation();

  // Read rigid transformation
  if (dofin_name != NULL) {
    transformation->irtkTransformation::Read(dofin_name);
  }

  // Set up registration and run
  registration->SetOutput(transformation);
  registration->Run();

  // Invert transformation
  if (invert == true) {
    transformation->Invert();
    transformation->UpdateParameter();
  }

  // Write rigid transformation
  if (dofout_name != NULL) {
    transformation->irtkTransformation::Write(dofout_name);
  }
}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] )
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
