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

irtkPointAffineRegistration::irtkPointAffineRegistration (): irtkPointRigidRegistration
    ()
{
  // Default optimization
  _OptimizationMethod = ConjugateGradientDescent;
}

irtkPointAffineRegistration::~irtkPointAffineRegistration ()
{
}

const char *irtkPointAffineRegistration::NameOfClass ()
{
  return "irtkPointAffineRegistration";
}

void irtkPointAffineRegistration::SetOutput (irtkTransformation * transformation)
{
  if (strcmp (transformation->NameOfClass (), "irtkAffineTransformation") != 0) {
    cerr << this->NameOfClass ()
         << "::SetOutput: Transformation must be affine" << endl;
    exit (0);
  }
  _transformation = transformation;
}

void irtkPointAffineRegistration::Run ()
{
  int i;
  irtkPoint p;
  irtkPointSet *tmp_target;
  irtkPointSet *tmp_source;
  irtkMatrix m(4, 4), m1(4, 4), m2(4, 4), m3(4, 4), m4(4, 4);

  // Do the initial set up
  this->Initialize();

  // Allocate memory for centered points
  irtkPointSet *centered_target = new irtkPointSet;
  irtkPointSet *centered_source = new irtkPointSet;

  // Calculate centroids
  irtkPoint centroid;
  for (i = 0; i < _target->Size (); i++) {
    p = _target->operator()(i);
    _transformation->Transform(p);
    centroid += p / _target->Size ();
  }

  // Calculate centred points
  for (i = 0; i < _target->Size (); i++) {
    p = _target->operator()(i);
    _transformation->Transform(p);
    centered_target->Add (p - centroid);
    p = _source->operator()(i);
    centered_source->Add (p - centroid);
  }

  // Get pointer to affine transformation
  irtkAffineTransformation *affine = (irtkAffineTransformation *)_transformation;

  // M1 is the current transformation matrix
  m1 = affine->GetMatrix();

  // Set affine transformation to identity
  m.Ident();
  affine->PutMatrix(m);
  affine->UpdateParameter();

  // Save old pointers
  tmp_target = _target;
  tmp_source = _source;

  // Setup new pointers
  _target = centered_target;
  _source = centered_source;

  // Optimize affine registration
  _optimizer->Run ();

  // Free memory
  delete centered_target;
  delete centered_source;

  // Restore old pointers
  _target = tmp_target;
  _source = tmp_source;

  // M2 is the centering transformation matrix
  m2.Ident();
  m2(0, 3) = -centroid._x;
  m2(1, 3) = -centroid._y;
  m2(2, 3) = -centroid._z;

  // M3 is the current affine transformation
  m3 = affine->GetMatrix();

  // M4 is the inverse centering matrix
  m4.Ident();
  m4(0, 3) = centroid._x;
  m4(1, 3) = centroid._y;
  m4(2, 3) = centroid._z;

  // Compute final transformation as M4 * M3 * M2 * M1
  affine->PutMatrix(m4 * m3 * m2 * m1);
  affine->UpdateParameter();

  // Do the final cleaning up
  this->Finalize ();

}
