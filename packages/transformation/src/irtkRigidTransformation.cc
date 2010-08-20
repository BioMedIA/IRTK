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

#ifdef HAS_VTK

extern Bool interactiveVTK;
extern Bool displayVTK;
extern Bool firstVTK;

extern void update   (irtkPointSet& pset, int nx, int ny, int nz);
extern void visualize(irtkPointSet& pset, int nx, int ny, int nz);

#endif

void irtkRigidTransformation::UpdateMatrix()
{
  // Update sines and cosines
  _cosrx = cos(_rx*(M_PI/180.0));
  _cosry = cos(_ry*(M_PI/180.0));
  _cosrz = cos(_rz*(M_PI/180.0));
  _sinrx = sin(_rx*(M_PI/180.0));
  _sinry = sin(_ry*(M_PI/180.0));
  _sinrz = sin(_rz*(M_PI/180.0));

  // Create a transformation whose transformation matrix is an identity matrix
  _matrix.Ident();

  // Add other transformation parameters to transformation matrix
  _matrix(0,0) = _cosry*_cosrz;
  _matrix(0,1) = _cosry*_sinrz;
  _matrix(0,2) = -_sinry;
  _matrix(0,3) = _tx;
  _matrix(1,0) = (_sinrx*_sinry*_cosrz-_cosrx*_sinrz);
  _matrix(1,1) = (_sinrx*_sinry*_sinrz+_cosrx*_cosrz);
  _matrix(1,2) = _sinrx*_cosry;
  _matrix(1,3) = _ty;
  _matrix(2,0) = (_cosrx*_sinry*_cosrz+_sinrx*_sinrz);
  _matrix(2,1) = (_cosrx*_sinry*_sinrz-_sinrx*_cosrz);
  _matrix(2,2) = _cosrx*_cosry;
  _matrix(2,3) = _tz;
  _matrix(3,3) = 1.0;
}

/// Construct a matrix based on parameters passed in the array.
irtkMatrix irtkRigidTransformation::Parameters2Matrix(double *params) const
{
  double tx = params[TX];
  double ty = params[TY];
  double tz = params[TZ];

  double rx = params[RX];
  double ry = params[RY];
  double rz = params[RZ];

  double cosrx = cos(rx*(M_PI/180.0));
  double cosry = cos(ry*(M_PI/180.0));
  double cosrz = cos(rz*(M_PI/180.0));
  double sinrx = sin(rx*(M_PI/180.0));
  double sinry = sin(ry*(M_PI/180.0));
  double sinrz = sin(rz*(M_PI/180.0));

  // Create a transformation whose transformation matrix is an identity matrix
  irtkMatrix matrix(4, 4);
  matrix.Ident();

  // Add other transformation parameters to transformation matrix
  matrix(0,0) = cosry*cosrz;
  matrix(0,1) = cosry*sinrz;
  matrix(0,2) = -sinry;
  matrix(0,3) = tx;

  matrix(1,0) = (sinrx*sinry*cosrz-cosrx*sinrz);
  matrix(1,1) = (sinrx*sinry*sinrz+cosrx*cosrz);
  matrix(1,2) = sinrx*cosry;
  matrix(1,3) = ty;

  matrix(2,0) = (cosrx*sinry*cosrz+sinrx*sinrz);
  matrix(2,1) = (cosrx*sinry*sinrz-sinrx*cosrz);
  matrix(2,2) = cosrx*cosry;
  matrix(2,3) = tz;
  matrix(3,3) = 1.0;

  return matrix;
}

void irtkRigidTransformation::UpdateParameter()
{
  double rigidParams[6];
  this->Matrix2Parameters(_matrix, rigidParams);

  _tx = rigidParams[TX];
  _ty = rigidParams[TY];
  _tz = rigidParams[TZ];

  _rx = rigidParams[RX];
  _ry = rigidParams[RY];
  _rz = rigidParams[RZ];

    // Update sines and cosines
  _cosrx = cos(_rx*(M_PI/180.0));
  _cosry = cos(_ry*(M_PI/180.0));
  _cosrz = cos(_rz*(M_PI/180.0));
  _sinrx = sin(_rx*(M_PI/180.0));
  _sinry = sin(_ry*(M_PI/180.0));
  _sinrz = sin(_rz*(M_PI/180.0));
}

void irtkRigidTransformation::Matrix2Parameters(irtkMatrix m, double* params) const
{
  double tmp;
  double TOL = 0.000001;

  params[TX] = m(0, 3);
  params[TY] = m(1, 3);
  params[TZ] = m(2, 3);

  tmp = asin(-1 * m(0, 2));

  // asin returns values for tmp in range -pi/2 to +pi/2, i.e. cos(tmp) >=
  // 0 so the division by cos(tmp) in the first part of the if clause was
  // not needed.
  if (fabs(cos(tmp)) > TOL) {
    params[RX] = atan2(m(1,2), m(2,2));
    params[RY] = tmp;
    params[RZ] = atan2(m(0,1), m(0,0));
  } else {
    //m(0,2) is close to +1 or -1
    params[RX] = atan2(-1.0*m(0,2)*m(1,0), -1.0*m(0,2)*m(2,0));
    params[RY] = tmp;
    params[RZ] = 0;
  }

  // Convert to degrees.
  params[RX] *= 180.0/M_PI;
  params[RY] *= 180.0/M_PI;
  params[RZ] *= 180.0/M_PI;

}

double irtkRigidTransformation::Get(int i) const
{
  switch (i) {
    case TX:
      return _tx;
      break;
    case TY:
      return _ty;
      break;
    case TZ:
      return _tz;
      break;
    case RX:
      return _rx;
      break;
    case RY:
      return _ry;
      break;
    case RZ:
      return _rz;
      break;
    default:
      cerr << "irtkRigidTransformation::Get: No such dof" << endl;
      return 0;
  }
}

void irtkRigidTransformation::Put(int i, double x)
{
  switch (i) {
    case TX:
      _tx = x;
      break;
    case TY:
      _ty = x;
      break;
    case TZ:
      _tz = x;
      break;
    case RX:
      _rx = x;
      break;
    case RY:
      _ry = x;
      break;
    case RZ:
      _rz = x;
      break;
    default:
      cerr << "irtkRigidTransformation::Put: No such dof" << endl;
      exit(1);
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

int irtkRigidTransformation::CheckHeader(char *name)
{
  int n;
  char buffer[255];

  ifstream from(name);

  if (!from) {
    cerr << "irtkRigidTransformation::CheckHeader: Can't open file "
    << name << "\n";
    exit(1);
  }

  // Read keyword
  from >> buffer;
  if (strcmp(buffer, "DOF:") != 0) {
    return False;
  }

  // Read no. of DOFs
  from >> n;
  if (n != 6) {
    return False;
  }

  return True;
}

void irtkRigidTransformation::Rotate(double& x, double& y, double& z)
{
  double a, b, c;

  // Pre-multiply point with transformation matrix
  a = _matrix(0, 0)*x+_matrix(0, 1)*y+_matrix(0, 2)*z;
  b = _matrix(1, 0)*x+_matrix(1, 1)*y+_matrix(1, 2)*z;
  c = _matrix(2, 0)*x+_matrix(2, 1)*y+_matrix(2, 2)*z;

  // Copy result back
  x = a;
  y = b;
  z = c;
}

void irtkRigidTransformation::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double)
{
  switch (dof) {
    case TX:
      jac[0]  = 1;
      jac[1]  = 0;
      jac[2]  = 0;
      break;
    case TY:
      jac[0]  = 0;
      jac[1]  = 1;
      jac[2]  = 0;
      break;
    case TZ:
      jac[0]  = 0;
      jac[1]  = 0;
      jac[2]  = 1;
      break;
    case RX:
      // Ensure that the derivatives are expressed with respect to angles not radians
      x *= (M_PI/180.0);
      y *= (M_PI/180.0);
      z *= (M_PI/180.0);
      jac[0]  = 0;
      jac[1]  = (_cosrx*_sinry*_cosrz+_sinrx*_sinrz)*x + (_cosrx*_sinry*_sinrz-_sinrx*_cosrz)*y + _cosrx*_cosry*z;
      jac[2]  = (-_sinrx*_sinry*_cosrz+_cosrx*_sinrz)*x + (-_sinrx*_sinry*_sinrz-_cosrx*_cosrz)*y - _sinrx*_cosry*z;
      break;
    case RY:
      // Ensure that the derivatives are expressed with respect to angles not radians
      x *= (M_PI/180.0);
      y *= (M_PI/180.0);
      z *= (M_PI/180.0);
      jac[0]  = -_sinry*_cosrz*x - _sinry*_sinrz*y - _cosry*z;
      jac[1]  = (_sinrx*_cosry*_cosrz)*x + (_sinrx*_cosry*_sinrz)*y - _sinry*_sinrx*z;
      jac[2]  = (_cosrx*_cosry*_cosrz)*x + (_cosrx*_cosry*_sinrz)*y - _cosrx*_sinry*z;
      break;
    case RZ:
      // Ensure that the derivatives are expressed with respect to angles not radians
      x *= (M_PI/180.0);
      y *= (M_PI/180.0);
      z *= (M_PI/180.0);
      jac[0]  = -_sinrz*_cosry*x + _cosry*_cosrz*y;
      jac[1]  = (-_sinrx*_sinry*_sinrz-_cosrx*_cosrz)*x + (_sinrx*_sinry*_cosrz-_cosrx*_sinrz)*y;
      jac[2]  = (-_cosrx*_sinry*_sinrz+_sinrx*_cosrz)*x + (_cosrx*_sinry*_cosrz+_sinrx*_sinrz)*y;
      break;
    default:
      cerr << this->NameOfClass() << "::JacobianDOFs(): No such dof = " << dof << endl;
      exit(1);
      break;
  }
}

void irtkRigidTransformation::Print()
{
  cout.setf(ios::right);
  cout.setf(ios::fixed);
  cout.precision(4);
  if (_status[TX]  == _Active) cout << "tx = " << setw(7) << _tx << " ";
  if (_status[TY]  == _Active) cout << "ty = " << setw(7) << _ty << " ";
  if (_status[TZ]  == _Active) cout << "tz = " << setw(7) << _tz << " ";
  if (_status[RX]  == _Active) cout << "rx = " << setw(7) << _rx << " ";
  if (_status[RY]  == _Active) cout << "ry = " << setw(7) << _ry << " ";
  if (_status[RZ]  == _Active) cout << "rz = " << setw(7) << _rz << endl;
  cout.precision(6);
  cout.unsetf(ios::right);
  cout.unsetf(ios::fixed);
}

Bool irtkRigidTransformation::IsIdentity()
{
  if ((_tx == 0) && (_ty == 0) && (_tz == 0) &&
      (_rx == 0) && (_ry == 0) && (_rz == 0)) {
    return True;
  } else {
    return False;
  }
}

istream& irtkRigidTransformation::Import(istream& is)
{
  int n;
  double dummy;
  char buffer[255];

  // Read keyword
  is >> buffer;
  if (strcmp(buffer, "DOF:") != 0) {
    cerr << "irtkRigidTransformation::Import: Not a valid transformation" << endl;
    exit(1);
  }

  // Read no. of DOFs
  is >> n;
  if (n != 6) {
    cerr << "irtkRigidTransformation::Import: Not a rigid transformation" << endl;
    exit(1);
  }

  // Read rigid transformation parameters
  is >> dummy >> dummy >> _tx;
  is >> dummy >> dummy >> _ty;
  is >> dummy >> dummy >> _tz;
  is >> dummy >> dummy >> _rx;
  is >> dummy >> dummy >> _ry;
  is >> dummy >> dummy >> _rz;

  // Update transformation matrix
  this->UpdateMatrix();

  // Invert transformation matrix and ...
  this->Invert();
  // ... and transformation parameters
  this->UpdateParameter();

  return is;
}

ostream& irtkRigidTransformation::Export(ostream& os)
{
  // Invert transformation matrix and ...
  this->Invert();
  // ... and transformation parameters
  this->UpdateParameter();

  // Write rigid transformation parameters
  os << "DOF: 6\n";
  os << "0.0\t0.0\t" << _tx << endl;
  os << "0.0\t0.0\t" << _ty << endl;
  os << "0.0\t0.0\t" << _tz << endl;
  os << "0.0\t0.0\t" << _rx << endl;
  os << "0.0\t0.0\t" << _ry << endl;
  os << "0.0\t0.0\t" << _rz << endl;

  // Invert transformation matrix and ...
  this->Invert();
  // ... and transformation parameters
  this->UpdateParameter();

  return os;
}

irtkCifstream& irtkRigidTransformation::Read(irtkCifstream& from)
{
  int i;
  double data;
  unsigned int magic_no, trans_type, dofs;

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkRigidTransformation::Read: Not a vaild transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if (trans_type != IRTKTRANSFORMATION_RIGID) {
    cerr << "irtkRigidTransformation::Read: Not a vaild rigid transformation " << trans_type << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&dofs, 1);
  if (dofs != 6) {
    cerr << "irtkRigidTransformation::Read: Invalid no. of dofs = " << dofs << endl;
    exit(1);
  }

  // Read data
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    from.ReadAsDouble(&data, 1);
    this->Put(i, data);
  }

  return from;
}

irtkCofstream& irtkRigidTransformation::Write(irtkCofstream& to)
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
