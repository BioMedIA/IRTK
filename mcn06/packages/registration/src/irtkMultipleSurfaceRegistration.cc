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

#include <irtkMultipleSurfaceRegistration.h>

#include <irtkUtil.h>

#define HISTORY

#ifdef HISTORY

extern irtkHistory *history;

#endif

irtkMultipleSurfaceRegistration::irtkMultipleSurfaceRegistration ()
{
  // Set inputs
  _target = NULL;
  _source = NULL;

  // Set output
  _transformation = NULL;

  // Others
  _target_locator = NULL;
  _source_locator = NULL;

  // Default parameters for registration
  _NumberOfIterations   = 100;
  _NumberOfSurfaces     = 0;
  _UseSymmetricDistance = False;
  _Epsilon              = 0.01;

#ifdef HISTORY
  history = new irtkHistory;
#endif

}

irtkMultipleSurfaceRegistration::~irtkMultipleSurfaceRegistration ()
{
#ifdef HISTORY
  delete history;
#endif
}

const char *irtkMultipleSurfaceRegistration::NameOfClass ()
{
  return "irtkMultipleSurfaceRegistration";
}

void irtkMultipleSurfaceRegistration::SetInput (vtkPolyData **target, vtkPolyData **source)
{
  _target = target;
  _source = source;
}

void irtkMultipleSurfaceRegistration::SetOutput (irtkTransformation *transformation)
{
  _transformation = transformation;
}

void irtkMultipleSurfaceRegistration::Initialize ()
{
  int i;

  if (_target == NULL) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Filter has no input for target" << endl;
    exit (1);
  }
  if (_source == NULL) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Filter has no input for source" << endl;
    exit (1);
  }
  if (_transformation == NULL) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Filter has no output" << endl;
    exit (1);
  }
  if (_source_locator == NULL) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Filter has no source locator" << endl;
    exit (1);
  }
  if (_NumberOfSurfaces == 0) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Filter has no surfaces" << endl;
    exit (1);
  }

  // Setup locator
  for (i = 0; i < _NumberOfSurfaces; i++) _source_locator[i]->SetDataSet(_source[i]);

  if ((_UseSymmetricDistance == True) && (_target_locator == NULL)) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Filter has no target locator" << endl;
    exit (1);
  } else {
    if (_UseSymmetricDistance == True) {
      // Setup locator
      for (i = 0; i < _NumberOfSurfaces; i++) _target_locator[i]->SetDataSet(_target[i]);
    }
  }
}

void irtkMultipleSurfaceRegistration::Finalize ()
{}

void irtkMultipleSurfaceRegistration::Optimize ()
{
  int i, j, k, n, id;
  double error, last_error, target_point[3], source_point[3], tmp_point[3];

  last_error = 0.0;

  for (j = 0; j < _NumberOfIterations; j++) {
    irtkPointSet source_pset, target_pset;
    n = 0;
    error = 0;
    for (k = 0; k < _NumberOfSurfaces; k++) {
      for (i = 0; i < _target[k]->GetNumberOfPoints(); i++) {
        _target[k]->GetPoints()->GetPoint (i, target_point);
        tmp_point[0] = target_point[0];
        tmp_point[1] = target_point[1];
        tmp_point[2] = target_point[2];
        _transformation->Transform(tmp_point[0], tmp_point[1], tmp_point[2]);
        source_point[0] = tmp_point[0];
        source_point[1] = tmp_point[1];
        source_point[2] = tmp_point[2];
        id = _source_locator[k]->FindClosestPoint (source_point);
        error += sqrt((tmp_point[0] - source_point[0]) * (tmp_point[0] - source_point[0]) +
                      (tmp_point[1] - source_point[1]) * (tmp_point[1] - source_point[1]) +
                      (tmp_point[2] - source_point[2]) * (tmp_point[2] - source_point[2]));
        target_pset.Add(target_point);
        source_pset.Add(source_point);
        n++;
      }
      if (_UseSymmetricDistance == True) {
        for (i = 0; i < _source[k]->GetNumberOfPoints(); i++) {
          _source[k]->GetPoints()->GetPoint (i, source_point);
          tmp_point[0] = source_point[0];
          tmp_point[1] = source_point[1];
          tmp_point[2] = source_point[2];
          _transformation->Inverse(tmp_point[0], tmp_point[1], tmp_point[2]);
          target_point[0] = tmp_point[0];
          target_point[1] = tmp_point[1];
          target_point[2] = tmp_point[2];
          id = _target_locator[k]->FindClosestPoint (target_point);
          error += sqrt((tmp_point[0] - target_point[0]) * (tmp_point[0] - target_point[0]) +
                        (tmp_point[1] - target_point[1]) * (tmp_point[1] - target_point[1]) +
                        (tmp_point[2] - target_point[2]) * (tmp_point[2] - target_point[2]));
          target_pset.Add(target_point);
          source_pset.Add(source_point);
          n++;
        }
      }
    }

    if (n == 0) {
      cout << "irtkMultipleSurfaceRegistration::Optimize: No corresponding points found" << endl;
      return;
    } else {
      error /= n;
    }

    cout << "Error = " << error << " at iteration = " << j << endl;

    if ((j > 0) && (last_error - error < _Epsilon)) {
      break;
    }
    last_error = error;

    // Run point registration
    _preg->SetInput(&target_pset, &source_pset);
    _preg->SetOutput(_transformation);
    _preg->Run();
  }
}

void irtkMultipleSurfaceRegistration::Run()
{
  // Initial setup
  this->Initialize();

  // Optimization
  this->Optimize();

  // Final setup
  this->Finalize();
}


#endif
