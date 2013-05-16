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

#include <irtkSurfaceRegistration.h>

#include <irtkUtil.h>

//#define HISTORY

irtkSurfaceRegistration::irtkSurfaceRegistration ()
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
  _ignore_edges         = false;
  _UseSymmetricDistance = false;
  _Epsilon              = 0.01;

#ifdef HISTORY
  history = new irtkHistory;
#endif

}

irtkSurfaceRegistration::~irtkSurfaceRegistration ()
{
#ifdef HISTORY
  delete history;
#endif
}

const char *irtkSurfaceRegistration::NameOfClass ()
{
  return "irtkSurfaceRegistration";
}

void irtkSurfaceRegistration::SetInput (vtkPolyData *target, vtkPolyData *source)
{
  _target = target;
  _source = source;
}

void irtkSurfaceRegistration::SetOutput (irtkTransformation *transformation)
{
  _transformation = transformation;
}

void irtkSurfaceRegistration::Initialize ()
{
  if (_source == NULL) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Filter has no input for source" << endl;
    exit (1);
  }
  if (_target == NULL) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Filter has no input for target" << endl;
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
  // Setup locator
  _source_locator->SetDataSet(_source);

  if ((_UseSymmetricDistance == true) && (_target_locator == NULL)) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Filter has no target locator" << endl;
    exit (1);
  } else {
    if (_UseSymmetricDistance == true) {
      // Setup locator
      _target_locator->SetDataSet(_target);
    }
  }
}

void irtkSurfaceRegistration::Finalize ()
{}

void irtkSurfaceRegistration::Optimize ()
{
  int i, j, n, id;
  double error, last_error, source_point[3], target_point[3], tmp_point[3];

  if (_ignore_edges && (_target->GetPointData ()->GetScalars () == NULL)) {
    cout << "irtkSurfaceRegistration::Optimize: Asked to ignore edges but no edge outline has been defined" << endl;
    exit (1);
  }

  last_error = 0.0;

  for (j = 0; j < _NumberOfIterations; j++) {
    irtkPointSet target_pset, source_pset;
    n = 0;
    error = 0;
    for (i = 0; i < _target->GetNumberOfPoints(); i++) {
      _target->GetPoints()->GetPoint (i, target_point);
      tmp_point[0] = target_point[0];
      tmp_point[1] = target_point[1];
      tmp_point[2] = target_point[2];
      _transformation->Transform(tmp_point[0], tmp_point[1], tmp_point[2]);
      source_point[0] = tmp_point[0];
      source_point[1] = tmp_point[1];
      source_point[2] = tmp_point[2];
      id = _source_locator->FindClosestPoint (source_point);
      if (_ignore_edges) {
        if (*_source->GetPointData()->GetScalars()->GetTuple (id)) {
          error += sqrt((tmp_point[0] - source_point[0]) * (tmp_point[0] - source_point[0]) +
                        (tmp_point[1] - source_point[1]) * (tmp_point[1] - source_point[1]) +
                        (tmp_point[2] - source_point[2]) * (tmp_point[2] - source_point[2]));
          target_pset.Add(target_point);
          source_pset.Add(source_point);
          n++;
        }
      } else {
        error += sqrt((tmp_point[0] - source_point[0]) * (tmp_point[0] - source_point[0]) +
                      (tmp_point[1] - source_point[1]) * (tmp_point[1] - source_point[1]) +
                      (tmp_point[2] - source_point[2]) * (tmp_point[2] - source_point[2]));
        target_pset.Add(target_point);
        source_pset.Add(source_point);
        n++;
      }
    }
    if (_UseSymmetricDistance == true) {
      for (i = 0; i < _source->GetNumberOfPoints(); i++) {
        _source->GetPoints()->GetPoint (i, source_point);
        tmp_point[0] = source_point[0];
        tmp_point[1] = source_point[1];
        tmp_point[2] = source_point[2];
        _transformation->Inverse(tmp_point[0], tmp_point[1], tmp_point[2]);
        target_point[0] = tmp_point[0];
        target_point[1] = tmp_point[1];
        target_point[2] = tmp_point[2];
        id = _target_locator->FindClosestPoint (target_point);
        if (_ignore_edges) {
          if (*_target->GetPointData()->GetScalars()->GetTuple (id)) {
            error += sqrt((tmp_point[0] - target_point[0]) * (tmp_point[0] - target_point[0]) +
                          (tmp_point[1] - target_point[1]) * (tmp_point[1] - target_point[1]) +
                          (tmp_point[2] - target_point[2]) * (tmp_point[2] - target_point[2]));
            target_pset.Add(target_point);
            source_pset.Add(source_point);
            n++;
          }
        } else {
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
      cout << "irtkSurfaceRegistration::Optimize: No corresponding points found" << endl;
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

void irtkSurfaceRegistration::Run()
{
  // Initial setup
  this->Initialize();

  // Optimization
  this->Optimize();

  // Final setup
  this->Finalize();
}


#endif
