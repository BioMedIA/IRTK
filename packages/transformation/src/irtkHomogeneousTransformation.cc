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

double irtkHomogeneousTransformation::Get(int index) const
{
  int i, j;

  if (index < this->NumberOfDOFs()) {
    i = index/4;
    j = index%4;
    return _matrix(i, j);
  } else {
    cerr << "irtkHomogeneousTransformation::Get: No such dof" << endl;
    return 0;
  }
}

void irtkHomogeneousTransformation::Put(int index, double x)
{
  int i, j;

  if (index < this->NumberOfDOFs()) {
    i = index/4;
    j = index%4;
    _matrix(i, j) = x;
  } else {
    cerr << "irtkHomogeneousTransformation::Put: No such dof" << endl;
    exit(1);
  }
}

void irtkHomogeneousTransformation::Import(char *name)
{
  // Open file
  ifstream from(name, ios::in | ios::binary);

  // Check whether file opened ok
  if (!from) {
    cerr << "irtkHomogeneousTransformation::Import: Can't open file "
         << name << "\n";
    exit(1);
  }

  this->Import(from);
}

void irtkHomogeneousTransformation::Export(char *name)
{
  // Open file
  ofstream to(name, ios::out | ios::binary);

  // Check whether file opened ok
  if (!to) {
    cerr << "irtkHomogeneousTransformation::Export: Can't open file "
         << name << "\n";
    exit(1);
  }

  this->Export(to);
}

void irtkHomogeneousTransformation::Print()
{
  _matrix.Print();
}

void irtkHomogeneousTransformation::Transform(double &x, double &y, double &z, double)
{
  double a, b, c;

  // Pre-multiply point with transformation matrix
  a = (_matrix(0, 0)*x+_matrix(0, 1)*y+_matrix(0, 2)*z+_matrix(0, 3));
  b = (_matrix(1, 0)*x+_matrix(1, 1)*y+_matrix(1, 2)*z+_matrix(1, 3));
  c = (_matrix(2, 0)*x+_matrix(2, 1)*y+_matrix(2, 2)*z+_matrix(2, 3));

  // Copy result back
  x = a;
  y = b;
  z = c;
}

void irtkHomogeneousTransformation::Displacement(double &x, double &y, double &z, double)
{
  double a, b, c;

  // Pre-multiply point with transformation matrix
  a = (_matrix(0, 0)*x+_matrix(0, 1)*y+_matrix(0, 2)*z+_matrix(0, 3));
  b = (_matrix(1, 0)*x+_matrix(1, 1)*y+_matrix(1, 2)*z+_matrix(1, 3));
  c = (_matrix(2, 0)*x+_matrix(2, 1)*y+_matrix(2, 2)*z+_matrix(2, 3));

  // Copy result back
  x = a - x;
  y = b - y;
  z = c - z;
}

double irtkHomogeneousTransformation::Inverse(double &x, double &y, double &z, double, double)
{
  double a, b, c;

  // Invert matrix
  irtkMatrix m = this->_matrix;
  m.Invert();

  // Pre-multiply point with transformation matrix
  a = (m(0, 0)*x+m(0, 1)*y+m(0, 2)*z+m(0, 3));
  b = (m(1, 0)*x+m(1, 1)*y+m(1, 2)*z+m(1, 3));
  c = (m(2, 0)*x+m(2, 1)*y+m(2, 2)*z+m(2, 3));

  // Copy result back
  x = a;
  y = b;
  z = c;

  // Error is zero by definition
  return 0;
}

Bool irtkHomogeneousTransformation::IsIdentity()
{
  int i, j;

  for (i = 0; i < _matrix.Rows(); i++) {
    for (j = 0; j < _matrix.Cols(); j++) {
      if (i == j) {
        if (_matrix(i, j) != 1) return False;
      } else {
        if (_matrix(i, j) != 0) return False;
      }
    }
  }
  return True;
}

void irtkHomogeneousTransformation::Invert()
{
  // Invert transformation
  _matrix.Invert();
}

void irtkHomogeneousTransformation::Jacobian(irtkMatrix &jac, double, double, double, double)
{
  int i, j;

  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Calculate matrix
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      jac(i, j) = _matrix(i, j);
    }
  }
}

void irtkHomogeneousTransformation::GlobalJacobian(irtkMatrix &jac, double, double, double, double)
{
  int i, j;

  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Calculate matrix
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      jac(i, j) = _matrix(i, j);
    }
  }
}

void irtkHomogeneousTransformation::LocalJacobian(irtkMatrix &jac, double, double, double, double)
{
  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Set matrix to identity
  jac(0, 0) = 1;
  jac(1, 1) = 1;
  jac(2, 2) = 1;
}


istream& irtkHomogeneousTransformation::Import(istream& is)
{
  is >> this->_matrix;

  // Invert transformation matrix
  this->Invert();

  return is;
}

ostream& irtkHomogeneousTransformation::Export(ostream& os)
{
  // Invert transformation matrix
  this->Invert();

  os << this->_matrix;

  // Invert transformation matrix
  this->Invert();

  return os;
}

irtkCifstream& irtkHomogeneousTransformation::Read(irtkCifstream& from)
{
  int i;
  double data;
  unsigned int magic_no, trans_type, dofs;

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkHomogeneousTransformation::Read: Not a vaild transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if (trans_type != IRTKTRANSFORMATION_RIGID) {
    cerr << "irtkHomogeneousTransformation::Read: Not a vaild affine transformation" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&dofs, 1);
  if (dofs != 12) {
    cerr << "irtkHomogeneousTransformation::Read: Invalid no. of dofs = " << dofs << endl;
    exit(1);
  }

  // Read data
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    from.ReadAsDouble(&data, 1);
    this->Put(i, data);
  }

  return from;
}

irtkCofstream& irtkHomogeneousTransformation::Write(irtkCofstream& to)
{
  int i;

  // Write magic no. for transformations
  unsigned int magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  unsigned int trans_type = IRTKTRANSFORMATION_RIGID;
  to.WriteAsUInt(&trans_type, 1);

  // Write transformation type
  unsigned int dofs = this->NumberOfDOFs();
  to.WriteAsUInt(&dofs, 1);

  // Write data
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    double data = this->Get(i);
    to.WriteAsDouble(&data, 1);
  }

  return to;
}
