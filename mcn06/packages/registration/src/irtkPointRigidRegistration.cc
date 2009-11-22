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

irtkPointRigidRegistration::irtkPointRigidRegistration (): irtkPointRegistration ()
{
  // Default optimization
  _OptimizationMethod = ClosedForm;
}

irtkPointRigidRegistration::~irtkPointRigidRegistration ()
{
}

char const *irtkPointRigidRegistration::NameOfClass ()
{
  return "irtkPointRigidRegistration";
}

void irtkPointRigidRegistration::SetOutput (irtkTransformation * transformation)
{
  if (strcmp (transformation->NameOfClass (), "irtkRigidTransformation") != 0) {
    cerr << this->NameOfClass ()
         << "::SetOutput: Transformation must be rigid" << endl;
    exit (0);
  }
  _transformation = transformation;
}

void irtkPointRigidRegistration::Initialize ()
{
  // Initialize base class
  this->irtkPointRegistration::Initialize ();

  // Initialize this class
  if (_source->Size () < 4) {
    cerr << this->NameOfClass ()
         << "::Initialize(): Must have at least four points" << endl;
    exit (1);
  }

}

void irtkPointRigidRegistration::ClosedFormOptimizer()
{
  int i;
  double x, y, z;
  irtkVector w;
  irtkMatrix u, v, r;

  // Local point sets
  irtkPointSet target (*_target);
  irtkPointSet source (*_source);

  // Calculate centroids
  irtkPoint target_centroid;
  irtkPoint source_centroid;
  for (i = 0; i < target.Size (); i++) {
    target_centroid -= target (i) / target.Size ();
    source_centroid -= source (i) / source.Size ();
  }

  // Subtract centroids
  for (i = 0; i < target.Size (); i++) {
    target (i) += target_centroid;
    source (i) += source_centroid;
  }

  // Calculate covariance matrix
  irtkMatrix h (3, 3);
  irtkMatrix a (3, 1);
  irtkMatrix b (1, 3);
  irtkMatrix m (4, 4);
  for (i = 0; i < target.Size (); i++) {
    a (0, 0) = target (i)._x;
    a (1, 0) = target (i)._y;
    a (2, 0) = target (i)._z;
    b (0, 0) = source (i)._x;
    b (0, 1) = source (i)._y;
    b (0, 2) = source (i)._z;
    h += a * b;
  }

  // Calculate SVD
  h.SVD (u, w, v);

  // Calculate rotation matrix
  u.Transpose ();
  r = v * u;

  // Check whether matrix is a rotation
  if (r.Det () > 0.999) {

    x = target_centroid._x;
    y = target_centroid._y;
    z = target_centroid._z;

    // Calculate rotated centroid
    target_centroid._x = r (0, 0) * x + r (0, 1) * y + r (0, 2) * z;
    target_centroid._y = r (1, 0) * x + r (1, 1) * y + r (1, 2) * z;
    target_centroid._z = r (2, 0) * x + r (2, 1) * y + r (2, 2) * z;

    // Calculate transformation matrix
    m (0, 0) = r (0, 0);
    m (1, 0) = r (1, 0);
    m (2, 0) = r (2, 0);
    m (3, 0) = 0;
    m (0, 1) = r (0, 1);
    m (1, 1) = r (1, 1);
    m (2, 1) = r (2, 1);
    m (3, 1) = 0;
    m (0, 2) = r (0, 2);
    m (1, 2) = r (1, 2);
    m (2, 2) = r (2, 2);
    m (3, 2) = 0;
    m (0, 3) = target_centroid._x - source_centroid._x;
    m (1, 3) = target_centroid._y - source_centroid._y;
    m (2, 3) = target_centroid._z - source_centroid._z;
    m (3, 3) = 1;

    // Update transformation
    ((irtkRigidTransformation *) _transformation)->PutMatrix (m);
    ((irtkRigidTransformation *) _transformation)->UpdateParameter();
  } else {

    cout << this->NameOfClass () << "::Run:Rotation involves reflection"
         << endl;

    // Search for most singular value
    i = 0;
    if ((w (0) < w (1)) && (w (0) < w (2)))
      i = 0;
    if ((w (1) < w (0)) && (w (1) < w (2)))
      i = 1;
    if ((w (2) < w (1)) && (w (2) < w (0)))
      i = 2;

    // Multiply column with most singular value by -1
    v (0, i) *= -1;
    v (1, i) *= -1;
    v (2, i) *= -1;

    // Recalculate rotation matrix
    r = v * u;

    x = target_centroid._x;
    y = target_centroid._y;
    z = target_centroid._z;

    // Calculate rotated centroid
    target_centroid._x = r (0, 0) * x + r (0, 1) * y + r (0, 2) * z;
    target_centroid._y = r (1, 0) * x + r (1, 1) * y + r (1, 2) * z;
    target_centroid._z = r (2, 0) * x + r (2, 1) * y + r (2, 2) * z;

    // Calculate transformation matrix
    m (0, 0) = r (0, 0);
    m (1, 0) = r (1, 0);
    m (2, 0) = r (2, 0);
    m (3, 0) = 0;
    m (0, 1) = r (0, 1);
    m (1, 1) = r (1, 1);
    m (2, 1) = r (2, 1);
    m (3, 1) = 0;
    m (0, 2) = r (0, 2);
    m (1, 2) = r (1, 2);
    m (2, 2) = r (2, 2);
    m (3, 2) = 0;
    m (0, 3) = target_centroid._x - source_centroid._x;
    m (1, 3) = target_centroid._y - source_centroid._y;
    m (2, 3) = target_centroid._z - source_centroid._z;
    m (3, 3) = 1;

    // Update transformation
    ((irtkRigidTransformation *) _transformation)->PutMatrix (m);
    ((irtkRigidTransformation *) _transformation)->UpdateParameter();
  }
}

void irtkPointRigidRegistration::Run ()
{
  // Do the initial set up
  this->Initialize ();

  // Optimize rigid registration
  if (_OptimizationMethod != ClosedForm) {
    cerr << "irtkPointRigidRegistration: Requires use of closed form optimization" << endl;
    exit(1);
  } else {
    this->ClosedFormOptimizer();
  }

  // Do the final cleaning up
  this->Finalize ();
}

