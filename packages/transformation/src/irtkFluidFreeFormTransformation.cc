/*=========================================================================
 
  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$
 
=========================================================================*/

#include <irtkTransformation.h>

#include <irtkFluidFreeFormTransformation.h>

#include <newt2.h>

irtkFluidFreeFormTransformation::irtkFluidFreeFormTransformation() : irtkMultiLevelFreeFormTransformation()
{}

irtkFluidFreeFormTransformation::irtkFluidFreeFormTransformation(const irtkFluidFreeFormTransformation &transformation) : irtkMultiLevelFreeFormTransformation(transformation)
{
  int i;

  cout << "Do not use" << endl;

  // Initialize local transformation
  for (i = 0; i < transformation._NumberOfLevels; i++) {
    _localTransformation[i] = transformation._localTransformation[i];
  }

  _NumberOfLevels = transformation._NumberOfLevels;
}

irtkFluidFreeFormTransformation::irtkFluidFreeFormTransformation(const irtkRigidTransformation &transformation) : irtkMultiLevelFreeFormTransformation(transformation)
{
  int i;

  // Initialize local transformation
  for (i = 0; i < MAX_TRANS; i++) {
    _localTransformation[i] = NULL;
  }

  _NumberOfLevels = 0;
}

irtkFluidFreeFormTransformation::irtkFluidFreeFormTransformation(const irtkAffineTransformation &transformation) : irtkMultiLevelFreeFormTransformation(transformation)
{
  int i;

  // Initialize local transformation
  for (i = 0; i < MAX_TRANS; i++) {
    _localTransformation[i] = NULL;
  }

  _NumberOfLevels = 0;
}


void irtkFluidFreeFormTransformation::Transform(double &x, double &y, double &z, double t)
{
  int i;

  // Compute local transformation
  for (i = _NumberOfLevels - 1; i >= 0; i--) {

    // Displacement
    _localTransformation[i]->Transform(x, y, z, t);
  }

  // Compute global transformation
  this->irtkAffineTransformation::Transform(x, y, z, t);

}

void irtkFluidFreeFormTransformation::Displacement(double& x, double& y, double& z, double t)
{
  double u, v, w;

  // Save current position
  u = x;
  v = y;
  w = z;

  // Transform current position
  this->Transform(u, v, w, t);

  x = u - x;
  y = v - y;
  z = w - z;
}

void irtkFluidFreeFormTransformation::Transform(int n, double &x, double &y, double &z, double t)
{
  int i;

  if (n > _NumberOfLevels) {
    cerr << "irtkFluidFreeFormTransformation::Transform: No such "
    << "transformation" << endl;
    exit(1);
  }

  // Compute local transformation
  for (i = n - 1; i >= 0; i--) {

    // Displacement
    _localTransformation[i]->Transform(x, y, z, t);
  }

  // Compute global transformation
  this->irtkAffineTransformation::Transform(x, y, z, t);
}

void irtkFluidFreeFormTransformation::GlobalTransform(double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::GlobalTransform: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::GlobalDisplacement(double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::GlobalDisplacement: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::LocalTransform(double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::LocalTransform: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::LocalTransform(int, double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::LocalTransform: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::LocalDisplacement(int, double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::LocalDisplacement: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::LocalDisplacement(double &, double &, double &, double)
{
  cerr << "irtkFluidFreeFormTransformation::LocalDisplacement: Does not make sense" << endl;
  exit(1);
}

void irtkFluidFreeFormTransformation::Jacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
  int i;
  irtkMatrix tmp_jac;

  // Initialize to identity
  jac.Ident();

  // Compute local jacobian
  for (i = _NumberOfLevels - 1; i >= 0; i--) {

    // Calculate jacobian
    _localTransformation[i]->Jacobian(tmp_jac, x, y, z, t);

    // Multiply jacobian
    jac = tmp_jac * jac;
  }

  // Compute global jacobian
  this->GlobalJacobian(tmp_jac, x, y, z, t);
  jac = tmp_jac * jac;
}

void irtkFluidFreeFormTransformation::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
  int i;
  irtkMatrix tmp_jac(3, 3);

  // Initialize to identity
  jac.Ident();

  // Compute local jacobian
  for (i = _NumberOfLevels - 1; i >= 0; i--) {

    // Calculate jacobian
    _localTransformation[i]->Jacobian(tmp_jac, x, y, z, t);

    // Multiply jacobian
    jac = tmp_jac * jac;
  }
}

irtkCifstream& irtkFluidFreeFormTransformation::Read(irtkCifstream& from)
{
  int i, offset;
  unsigned int magic_no, trans_type;

  // Delete old local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    if (_localTransformation[i] != NULL) delete _localTransformation[i];
  }

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkFluidFreeFormTransformation::Read: Not a vaild transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if (trans_type != IRTKTRANSFORMATION_FLUID) {
    cerr << "irtkFluidFreeFormTransformation::Read: Not a vaild fluid transformation" << endl;
    exit(1);
  }

  // Read number of local transformations
  from.ReadAsInt(&_NumberOfLevels, 1);

  // Read global transformation
  this->irtkAffineTransformation::Read(from);

  // Read local transformations
  for (i = 0; i < _NumberOfLevels; i++) {

    // Remember current file location
    offset = from.Tell();

    // Read magic no. for transformations
    from.ReadAsUInt(&magic_no, 1);
    if (magic_no != IRTKTRANSFORMATION_MAGIC) {
      cerr << "irtkFluidFreeFormTransformation::Read: Not a valid fluid transformation found at = " << offset << endl;
      exit(1);
    }

    // Read transformation type
    from.ReadAsUInt(&trans_type, 1);
    if ((trans_type == IRTKTRANSFORMATION_LINEAR_FFD) || (trans_type == IRTKTRANSFORMATION_LINEAR_FFD_EXT1)) {
      from.Seek(offset);
      _localTransformation[i] = new irtkLinearFreeFormTransformation;
      ((irtkLinearFreeFormTransformation *)_localTransformation[i])->Read(from);
    } else {
      if ((trans_type == IRTKTRANSFORMATION_BSPLINE_FFD) || (trans_type == IRTKTRANSFORMATION_BSPLINE_FFD_EXT1)) {
        from.Seek(offset);
        _localTransformation[i] = new irtkBSplineFreeFormTransformation;
        ((irtkBSplineFreeFormTransformation *)_localTransformation[i])->Read(from);
      } else {
        cerr << "irtkFluidFreeFormTransformation::Read: No a valid FFD transformation type at = " << offset << " " << endl;
        exit(1);
      }
    }
  }
  return from;
}

irtkCofstream& irtkFluidFreeFormTransformation::Write(irtkCofstream& to)
{
  int i;
  unsigned int magic_no, trans_type;

  // Write magic no. for transformations
  magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  trans_type = IRTKTRANSFORMATION_FLUID;
  to.WriteAsUInt(&trans_type, 1);

  // Write transformation type
  to.WriteAsInt(&_NumberOfLevels, 1);

  // Write global transformation
  this->irtkAffineTransformation::Write(to);

  // Write local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    _localTransformation[i]->Write(to);
  }

  return to;
}

double irtkFluidFreeFormTransformation::Inverse(double &x, double &y, double &z, double t, double tolerance)
{
  int i;
  double error;

  // Invert global transformation
  error = this->irtkHomogeneousTransformation::Inverse(x, y, z, t);

  // Invert local transformations
  for (i = 0; i < _NumberOfLevels; i++) {
    error += _localTransformation[i]->Inverse(x, y, z, t, tolerance);
  }

  // Return error
  return error;
}

double irtkFluidFreeFormTransformation::Inverse(int n, double &x, double &y, double &z, double t, double tolerance)
{
  int i;
  double error;

  if (n > _NumberOfLevels) {
    cerr << "irtkFluidFreeFormTransformation::Inverse: No such "
    << "transformation" << endl;
    exit(1);
  }

  // Invert global transformation
  error = this->irtkHomogeneousTransformation::Inverse(x, y, z, t);

  // Invert local transformations
  for (i = 0; i < n; i++) {
    error += _localTransformation[i]->Inverse(x, y, z, t, tolerance);
  }

  // Return error
  return error;
}
