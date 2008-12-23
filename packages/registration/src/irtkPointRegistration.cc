/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

#ifdef HAS_VTK
#include <vtkPolyDataWriter.h>
#endif

irtkPointRegistration::irtkPointRegistration ()
{
  // Set inputs
  _target = NULL;
  _source = NULL;

  // Set output
  _transformation = NULL;
}

irtkPointRegistration::~irtkPointRegistration ()
{
}

const char *irtkPointRegistration::NameOfClass ()
{
  return "irtkPointRegistration";
}

void irtkPointRegistration::SetInput (irtkPointSet * target, irtkPointSet * source)
{
  _target = target;
  _source = source;
}

void irtkPointRegistration::SetOutput (irtkTransformation * transformation)
{
  _transformation = transformation;
}

void irtkPointRegistration::Initialize ()
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

  // Setup the optimizer
  switch (_OptimizationMethod) {
  case DownhillDescent:
    _optimizer = new irtkDownhillDescentOptimizer;
    break;
  case GradientDescent:
    _optimizer = new irtkGradientDescentOptimizer;
    break;
  case SteepestGradientDescent:
    _optimizer = new irtkSteepestGradientDescentOptimizer;
    break;
  case ConjugateGradientDescent:
    _optimizer = new irtkConjugateGradientDescentOptimizer;
    break;
  case ClosedForm:
    _optimizer = NULL;
    break;
  default:
    cerr << "Unkown optimizer" << endl;
    exit (1);
  }

  if (_optimizer != NULL) {
    _optimizer->SetTransformation (_transformation);
    _optimizer->SetRegistration (this);

  }
}

void irtkPointRegistration::Finalize ()
{
  if (_optimizer != NULL) delete _optimizer;
}

double irtkPointRegistration::Evaluate ()
{
  int i;
  double distance = 0;
  irtkPoint _temp_b;

  for (i = 0; i < _target->Size (); i++) {
    irtkPoint a = _target->operator()(i);
    irtkPoint b = _source->operator()(i);
    _transformation->Transform(a);
    distance += a.Distance (b);
  }
  return -distance;
}

double irtkPointRegistration::EvaluateGradient (float step, float *dx)
{
  int i;
  double s1, s2, norm, parameterValue;

  for (i = 0; i < _transformation->NumberOfDOFs (); i++) {
    if (_transformation->irtkTransformation::GetStatus (i) == _Active) {
      parameterValue = _transformation->Get (i);
      _transformation->Put (i, parameterValue + step);
      s1 = this->Evaluate ();
      _transformation->Put (i, parameterValue - step);
      s2 = this->Evaluate ();
      _transformation->Put (i, parameterValue);
      dx[i] = s1 - s2;
    } else {
      dx[i] = 0;
    }
  }

  // Calculate norm of vector
  norm = 0;
  for (i = 0; i < _transformation->NumberOfDOFs (); i++) {
    norm += dx[i] * dx[i];
  }

  // Normalize vector
  norm = sqrt (norm);
  if (norm > 0) {
    for (i = 0; i < _transformation->NumberOfDOFs (); i++) {
      dx[i] /= norm;
    }
  } else {
    for (i = 0; i < _transformation->NumberOfDOFs (); i++) {
      dx[i] = 0;
    }
  }
  return norm;
}
