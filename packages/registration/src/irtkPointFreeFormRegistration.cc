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

irtkPointFreeFormRegistration::irtkPointFreeFormRegistration (): irtkPointRegistration
    ()
{
}

irtkPointFreeFormRegistration::~irtkPointFreeFormRegistration ()
{
}

const char *irtkPointFreeFormRegistration::NameOfClass ()
{
  return "irtkPointFreeFormRegistration";
}

void irtkPointFreeFormRegistration::SetOutput( irtkTransformation *transformation)
{
  if (strcmp (transformation->NameOfClass (), "irtkMultiLevelFreeFormTransformation") != 0) {
    cerr << this->NameOfClass ()
         << "::SetOutput: Transformation must be irtkMultiLevelFreeFormTransformation" << endl;
    exit (0);
  }
  _transformation = transformation;
}

void irtkPointFreeFormRegistration::Initialize ()
{
  // Initialize this class
  if (_source != NULL) {
    if (_source->Size () < 4) {
      cerr << this->NameOfClass ()
           << "::Initialize(): Must have at least four points" << endl;
      exit (1);
    }
  }

  // Pointer to multi-level FFD
  _mffd = (irtkMultiLevelFreeFormTransformation *)_transformation;

  // Pointer to single-level FFD
  _affd = dynamic_cast<irtkBSplineFreeFormTransformation *>(_mffd->PopLocalTransformation());
  if (_affd == NULL) {
    cerr << "irtkPointFreeFormRegistration::Initialize: Current level is not a B-spline transformation" << endl;
    exit(1);
  }
}

void irtkPointFreeFormRegistration::Finalize()
{
  // Push local transformation back on transformation stack
  _mffd->PushLocalTransformation(_affd);
}

void irtkPointFreeFormRegistration::Optimize ()
{
  double *x_target;
  double *y_target;
  double *z_target;

  double *x_source;
  double *y_source;
  double *z_source;

  irtkPoint target_point, source_point;

  if (_target != NULL) {
    x_target = new double[_target->Size ()];
    y_target = new double[_target->Size ()];
    z_target = new double[_target->Size ()];

    x_source = new double[_source->Size ()];
    y_source = new double[_source->Size ()];
    z_source = new double[_source->Size ()];

    for (int count (0); count < _source->Size (); count++) {
      target_point = (*_target) (count);
      x_target[count] = target_point._x;
      y_target[count] = target_point._y;
      z_target[count] = target_point._z;

      // Transform points
      _mffd->Transform(target_point._x, target_point._y, target_point._z);

      // Calculate residual distance
      source_point = (*_source) (count);
      x_source[count] = source_point._x - target_point._x;
      y_source[count] = source_point._y - target_point._y;
      z_source[count] = source_point._z - target_point._z;
    }

    // Approximate residual difference
    _affd->Approximate (x_target, y_target, z_target, x_source,
                        y_source, z_source, _target->Size ());

    delete []x_target;
    delete []y_target;
    delete []z_target;
    delete []x_source;
    delete []y_source;
    delete []z_source;
  }
}

void irtkPointFreeFormRegistration::Run ()
{
  // Do the initial set up
  this->Initialize ();

  // Optimize rigid registration
  this->Optimize ();

  // Do the final cleaning up
  this->Finalize ();
}

