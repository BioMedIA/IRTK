/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#ifdef HAS_VTK

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataWriter.h>

#include <irtkTransformation.h>

#include <irtkSurfaceRegistration.h>

char *_target_name = NULL, *_source_name = NULL;
char *dofin_name = NULL, *dofout_name = NULL;

void usage()
{
  cerr << "Usage: sareg [target] [source] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-locator>          Locator: 0 = cell locator, 1 = point locator, 2 = kd-tree locator (default = 1)" << endl;
  cerr << "<-dofin name>       Name of input file" << endl;
  cerr << "<-dofout name>      Name of output file" << endl;
  cerr << "<-epsilon>          Value for espilon (default=0.01)" << endl;
  cerr << "<-clean>            Clean polydata (default OFF)" << endl;
  cerr << "<-symmetric>        Use symmetric distance (default OFF)" << endl;
  cerr << "<-ignoreedges>      Ignores edges in ICP (default OFF)" << endl;
  cerr << "<-invert>           Save the inverse transformation (default OFF)" << endl;
  cerr << "<-iterations>       Number of 3D registration iterations (default 100)" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  float epsilon;
  int locatorType, iterations, ok;
  bool clean, invert, ignoreEdges, symmetricDistance;

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
  ignoreEdges = false;
  symmetricDistance = false;

  // Parse filenames
  _target_name = argv[1];
  argv++;
  argc--;
  _source_name = argv[1];
  argv++;
  argc--;

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
    if ((ok == false) && (strcmp(argv[1], "-ignoreedges") == 0)) {
      argc--;
      argv++;
      ignoreEdges = true;
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
    if ((ok == false) && (strcmp(argv[1], "-clean") == 0)) {
      argc--;
      argv++;
      clean = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-symmetric") == 0)) {
      argc--;
      argv++;
      symmetricDistance = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
      argc--;
      argv++;
      invert = true;
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

  // Target pipeline
  cout << "Reading target ... " << _target_name << endl;
  vtkPolyDataReader *target_reader = vtkPolyDataReader::New();
  target_reader->SetFileName(_target_name);
  target_reader->Modified();
  target_reader->Update();
  vtkPolyData *target = vtkPolyData::New();

  if (clean == true) {
    vtkCleanPolyData *target_cleaner = vtkCleanPolyData::New();
#if VTK_MAJOR_VERSION >= 6
        target_cleaner->SetInputData(target_reader->GetOutput());
#else //VTK_MAJOR_VERSION >= 6
        target_cleaner->SetInput(target_reader->GetOutput());
#endif //VTK_MAJOR_VERSION >= 6
    target_cleaner->Modified();
    target_cleaner->Update();
    target = target_cleaner->GetOutput();
#if VTK_MAJOR_VERSION < 6
    target->Update();
#endif

  } else {
    target = target_reader->GetOutput();
#if VTK_MAJOR_VERSION < 6
    target->Update();
#endif
  }

  // Source pipeline
  cout << "Reading source ... " << _source_name << endl;
  vtkPolyDataReader *source_reader = vtkPolyDataReader::New();
  source_reader->SetFileName(_source_name);
  source_reader->Modified();
  source_reader->Update();
  vtkPolyData *source = vtkPolyData::New();

  if (clean == true) {
    vtkCleanPolyData *source_cleaner = vtkCleanPolyData::New();
#if VTK_MAJOR_VERSION >= 6
    source_cleaner->SetInputData(source_reader->GetOutput());
#else //VTK_MAJOR_VERSION >= 6
    source_cleaner->SetInput(source_reader->GetOutput());
#endif //VTK_MAJOR_VERSION >= 6

    source_cleaner->Modified();
    source_cleaner->Update();
    source = source_cleaner->GetOutput();
#if VTK_MAJOR_VERSION < 6
    source->Update();
#endif

  } else {
    source = source_reader->GetOutput();
#if VTK_MAJOR_VERSION < 6
    source->Update();
#endif
  }

  if (ignoreEdges == true) {
    MarkBoundary(target);
    MarkBoundary(source);
  }

  // Create registration
  irtkSurfaceAffineRegistration *registration = new  irtkSurfaceAffineRegistration;
  registration->SetInput(target, source);
  registration->SetNumberOfIterations(iterations);
  registration->SetEpsilon(epsilon);

  // Check if to do symmetric registration
  if (symmetricDistance == true) {
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

  // Check if to ignore edges
  if (ignoreEdges == true) {
    registration->IgnoreEdges();
  }

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
