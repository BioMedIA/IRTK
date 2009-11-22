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

void irtkAffineTransformation::UpdateMatrix()
{
  // Update rigid transformation
  irtkRigidTransformation::UpdateMatrix();

  // Update affine transformation: Add shearing
  irtkMatrix skew(4, 4);
  skew.Ident();
  skew(0, 1) = tan(_sxy*(M_PI/180.0));
  skew(0, 2) = tan(_sxz*(M_PI/180.0));
  skew(1, 2) = tan(_syz*(M_PI/180.0));
  _matrix *= skew;

  // Update affine transformation: Add scaling
  irtkMatrix scale(4, 4);
  scale.Ident();
  scale(0, 0) = _sx / 100.0;
  scale(1, 1) = _sy / 100.0;
  scale(2, 2) = _sz / 100.0;
  _matrix *= scale;

}

/// Construct a matrix based on parameters passed in the array.
irtkMatrix irtkAffineTransformation::Parameters2Matrix(double *params) const
{
  // Get the rigid components.
  irtkMatrix matrix;
  matrix = this->irtkRigidTransformation::Parameters2Matrix(params);

  // Read scale and shear parameters.
  double sx  = params[SX];
  double sy  = params[SY];
  double sz  = params[SZ];

  double sxy = params[SXY];
  double syz = params[SYZ];
  double sxz = params[SXZ];

  // Update affine transformation: Add shearing
  irtkMatrix skew(4, 4);
  skew.Ident();
  skew(0, 1) = tan(sxy*(M_PI/180.0));
  skew(0, 2) = tan(sxz*(M_PI/180.0));
  skew(1, 2) = tan(syz*(M_PI/180.0));
  matrix *= skew;

  // Update affine transformation: Add scaling
  irtkMatrix scale(4, 4);
  scale.Ident();
  scale(0, 0) = sx / 100.0;
  scale(1, 1) = sy / 100.0;
  scale(2, 2) = sz / 100.0;
  matrix *= scale;

  return matrix;
}

/// Construct a matrix based on old style parameters.
irtkMatrix irtkAffineTransformation::Parameters2Matrix_15DOFs(double *params) const
{
  // Get the rigid components.
  irtkMatrix matrix;
  matrix = this->irtkRigidTransformation::Parameters2Matrix(params);

  // Read scale and shear parameters.
  double sx  = params[SX];
  double sy  = params[SY];
  double sz  = params[SZ];

  double sxy = params[SXY];
  double syx = params[SYX];
  double syz = params[SYZ];
  double szy = params[SZY];
  double szx = params[SZX];
  double sxz = params[SXZ];

  // Update affine transformation: Add scaling
  irtkMatrix scale(4, 4);
  scale.Ident();
  scale(0, 0) = sx / 100.0;
  scale(1, 1) = sy / 100.0;
  scale(2, 2) = sz / 100.0;
  matrix *= scale;

  // Update affine transformation: Add shearing
  irtkMatrix skewx(4, 4);
  skewx.Ident();
  skewx(2, 1)= tan(szy*(M_PI/180.0));
  skewx(1, 2)= tan(syz*(M_PI/180.0));
  irtkMatrix skewy(4, 4);
  skewy.Ident();
  skewy(2, 0)= tan(szx*(M_PI/180.0));
  skewy(0, 2)= tan(sxz*(M_PI/180.0));
  irtkMatrix skewz(4, 4);
  skewz.Ident();
  skewz(1, 0)= tan(sxy*(M_PI/180.0));
  skewz(0, 1)= tan(syx*(M_PI/180.0));
  matrix *= skewx * skewy * skewz;

  return matrix;
}

// Use current matrix to extract the transformation parameters based on 12
// DOF model, see:
//    http://www.acm.org/pubs/tog/GraphicsGems/gemsii/unmatrix.c
// (from Graphics Gems II)
void irtkAffineTransformation::UpdateParameter()
{
  // Use _matrix to evaluate the affine transformation parameters.
  // It is assumed that there is no perspective transformation.
  int i;
  double TOL = 0.000001;
  double tansxy, tansxz, tansyz;

  if (fabs(_matrix(3,3) - 1.0) > TOL) {
    cerr << "irtkAffineTransformation::UpdateParameter" << endl;
    cerr << "Value at _matrix(3,3) must equal 1." << endl;
    exit(1);
  }

  if (fabs(_matrix.Det()) < TOL) {
    cerr << "irtkAffineTransformation::UpdateParameter" << endl;
    cerr << "Matrix singular (or very close to singular!)." << endl;
    exit(1);
  }

  // First Part Of Graphics Gems Code Ignored Because It Relates To
  // Perspective Transformation.
  if (fabs(_matrix(3, 0)) > TOL ||
      fabs(_matrix(3, 1)) > TOL ||
      fabs(_matrix(3, 2)) > TOL ) {
    cerr << "irtkAffineTransformation::UpdateParameter" << endl;
    cerr << "Matrix contains perspective components." << endl;
    exit(1);
  }

  irtkMatrix copy(4, 4);
  copy = _matrix;

  // Get scale and shear by manipulating the columns of the upper left 3x3
  // sub-matrix.
  irtkVector col_0, col_1, col_2;
  col_0.Initialize(3);
  col_1.Initialize(3);
  col_2.Initialize(3);
  for (i = 0; i < 3; ++i) {
    col_0(i) = copy(i, 0);
    col_1(i) = copy(i, 1);
    col_2(i) = copy(i, 2);
  }

  // Compute X scale factor and normalize first col.
  _sx = col_0.Norm();
  col_0 /= _sx;

  // Compute XY shear factor and make 2nd col orthogonal to 1st.
  tansxy = col_0.ScalarProduct(col_1);
  col_1 = col_1 - col_0 * tansxy;

  // Actually, tansxy and col_1 are still to large by a factor of sy.
  // Now, compute Y scale and normalize 2nd col and rescale tansxy.
  _sy = col_1.Norm();
  col_1  /= _sy;
  tansxy /= _sy;

  // Compute XZ and YZ shears, orthogonalize 3rd col
  tansxz = col_0.ScalarProduct(col_2);
  col_2 = col_2 - col_0 * tansxz;

  tansyz = col_1.ScalarProduct(col_2);
  col_2 = col_2 - col_1 * tansyz;

  // Actually, tansxz, tansyz and col_2 are still to large by a factor of
  // sz.  Next, get Z scale, normalize 3rd col and scale tansxz and tansyz.
  _sz = col_2.Norm();
  col_2  /= _sz;
  tansxz /= _sz;
  tansyz /= _sz;

  // At this point, the columns are orthonormal.  Check for a coordinate
  // system flip.  If the determinant is -1, then negate the matrix and the
  // scaling factors.
  irtkVector col_1_x_col_2;
  col_1_x_col_2.Initialize(3);
  col_1_x_col_2 = col_1.CrossProduct(col_2);

  if (col_0.ScalarProduct(col_1_x_col_2) < 0) {
    _sx *= -1;
    _sy *= -1;
    _sz *= -1;
    col_0 *= -1;
    col_1 *= -1;
    col_2 *= -1;
  }

  // Retrieve the shear angles in degrees.
  _sxy = atan(tansxy) * 180.0 / M_PI ;
  _sxz = atan(tansxz) * 180.0 / M_PI ;
  _syz = atan(tansyz) * 180.0 / M_PI ;

  // Convert scales to percentages.
  _sx *= 100;
  _sy *= 100;
  _sz *= 100;

  // Now get the rigid transformation parameters.
  double rigidParams[6];
  // Put the rotation matrix components into the upper left 3x3 submatrix.
  for (i = 0; i < 3; ++i) {
    copy(i, 0) = col_0(i);
    copy(i, 1) = col_1(i);
    copy(i, 2) = col_2(i);
  }

  this->irtkRigidTransformation::Matrix2Parameters(copy, rigidParams);

  _tx = rigidParams[TX];
  _ty = rigidParams[TY];
  _tz = rigidParams[TZ];
  _rx = rigidParams[RX];
  _ry = rigidParams[RY];
  _rz = rigidParams[RZ];
}

double irtkAffineTransformation::Get(int i) const
{
  if (i >= this->NumberOfDOFs()) {
    cerr << "irtkAffineTransformation::Get: No such dof" << endl;
    return 0;
  }
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
  case SX:
    return _sx;
    break;
  case SY:
    return _sy;
    break;
  case SZ:
    return _sz;
    break;
  case SXY:
    return _sxy;
    break;
  case SYZ:
    return _syz;
    break;
  case SXZ:
    return _sxz;
    break;
  default:
    cerr << "irtkAffineTransformation::Get: No such dof" << endl;
    return 0;
  }
}

void irtkAffineTransformation::Put(int i, double x)
{
  if (i >= this->NumberOfDOFs()) {
    cerr << "irtkAffineTransformation::Put: No such dof" << endl;
    return;
  }
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
  case SX:
    _sx = x;
    break;
  case SY:
    _sy = x;
    break;
  case SZ:
    _sz = x;
    break;
  case SXY:
    _sxy = x;
    break;
  case SYZ:
    _syz = x;
    break;
  case SXZ:
    _sxz = x;
    break;
  default:
    cerr << "irtkAffineTransformation::Put(): No such dof" << endl;
    exit(1);
  }
  // Update transformation matrix
  this->UpdateMatrix();
}

int irtkAffineTransformation::CheckHeader(char *name)
{
  int n;
  char buffer[255];

  ifstream from(name);

  if (!from) {
    cerr << "irtkAffineTransformation::CheckHeader: Can't open file "
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
  if ((n != 6) && (n != 9) && (n != 12) && (n != 15)) {
    return False;
  }

  return True;
}

void irtkAffineTransformation::Print()
{
  cout.setf(ios::right);
  cout.setf(ios::fixed);
  cout.precision(4);
  if (_status[TX]  == _Active) cout << "tx = "  << setw(7) << _tx << " ";
  if (_status[TY]  == _Active) cout << "ty = "  << setw(7) << _ty << " ";
  if (_status[TZ]  == _Active) cout << "tz = "  << setw(7) << _tz << " ";
  if (_status[RX]  == _Active) cout << "rx = "  << setw(7) << _rx << " ";
  if (_status[RY]  == _Active) cout << "ry = "  << setw(7) << _ry << " ";
  if (_status[RZ]  == _Active) cout << "rz = "  << setw(7) << _rz << " ";
  cout << endl;
  if (_status[SX]  == _Active) cout << "sx = "  << setw(7) << _sx << " ";
  if (_status[SY]  == _Active) cout << "sy = "  << setw(7) << _sy << " ";
  if (_status[SZ]  == _Active) cout << "sz = "  << setw(7) << _sz << " ";
  if (_status[SXY] == _Active) cout << "sxy = " << setw(7) << _sxy << " ";
  if (_status[SYZ] == _Active) cout << "syz = " << setw(7) << _syz << " ";
  if (_status[SXZ] == _Active) cout << "sxz = " << setw(7) << _sxz << " ";
  cout << endl;

  cout.precision(6);
  cout.unsetf(ios::right);
  cout.unsetf(ios::fixed);
}

Bool irtkAffineTransformation::IsIdentity()
{
  if ((_tx == 0)  && (_ty == 0)  && (_tz == 0)  &&
      (_rx == 0)  && (_ry == 0)  && (_rz == 0)  &&
      (_sx == 100.0)  && (_sy == 100.0)  && (_sz == 100.0)  &&
      (_sxy == 0) && (_syz == 0) && (_sxz == 0)) {
    return True;
  } else {
    return False;
  }
}

istream& irtkAffineTransformation::Import(istream& is)
{
  int n;
  double dummy;
  char buffer[255];
  double params[15];
  memset ((void*)&params[0], 0, 15 * sizeof(double));

  // Read keyword
  is >> buffer;
  if (strcmp(buffer, "DOF:") != 0) {
    cerr << "irtkAffineTransformation::operator>>: Not a valid transformation"
         << endl;
    exit(1);
  }

  // Read no. of DOFs
  is >> n;
  if ((n != 6)  && (n != 9) && (n != 12) && (n != 15)) {
    cerr << "irtkAffineTransformation::operator>>: Not an affine or rigid "
         << "transformation" << endl;
    exit(1);
  }

  // Read rigid transformation parameters
  is >> dummy >> dummy >> params[TX];
  is >> dummy >> dummy >> params[TY];
  is >> dummy >> dummy >> params[TZ];
  is >> dummy >> dummy >> params[RX];
  is >> dummy >> dummy >> params[RY];
  is >> dummy >> dummy >> params[RZ];

  if (n > 6) {
    // Read affine transformation parameters (scaling)
    is >> dummy >> dummy >> params[SX]; params[SX] *= 100;
    is >> dummy >> dummy >> params[SY]; params[SY] *= 100;
    is >> dummy >> dummy >> params[SZ]; params[SZ] *= 100;
  } else {
    params[SX] = 100;
    params[SY] = 100;
    params[SZ] = 100;
  }

  if (n == 12) {
    is >> dummy >> dummy >> params[SXY];
    is >> dummy >> dummy >> params[SYZ];
    is >> dummy >> dummy >> params[SXZ];
  } else if (n == 15) {
    // Old style, 6 skew parameters.
    is >> dummy >> dummy >> params[SXY];
    is >> dummy >> dummy >> params[SYX];
    is >> dummy >> dummy >> params[SYZ];
    is >> dummy >> dummy >> params[SZY];
    is >> dummy >> dummy >> params[SZX];
    is >> dummy >> dummy >> params[SXZ];
  }

  if (n <= 12) {
    // Call to set params also updates transformation matrix.
    this->SetParameters(&params[0]);
  } else { // n == 15
    irtkMatrix temp = this->Parameters2Matrix_15DOFs(&params[0]);
    // Call to put matrix also updates parameters in a 12DOF format.
    this->PutMatrix(temp);
  }

  // Invert transformation matrix and ...
  this->Invert();
  // ... and transformation parameters
  this->UpdateParameter();

  return is;
}

ostream& irtkAffineTransformation::Export(ostream& os)
{
  // Invert transformation matrix and ...
  this->Invert();
  // ... and transformation parameters
  this->UpdateParameter();

  // Write affine transformation parameters
  os << "DOF: " << this->NumberOfDOFs() << endl;
  os << "0.0\t0.0\t" << _tx << endl;
  os << "0.0\t0.0\t" << _ty << endl;
  os << "0.0\t0.0\t" << _tz << endl;
  os << "0.0\t0.0\t" << _rx << endl;
  os << "0.0\t0.0\t" << _ry << endl;
  os << "0.0\t0.0\t" << _rz << endl;
  os << "0.0\t0.0\t" << _sx/100.0 << endl;
  os << "0.0\t0.0\t" << _sy/100.0 << endl;
  os << "0.0\t0.0\t" << _sz/100.0 << endl;
  os << "0.0\t0.0\t" << _sxy << endl;
  os << "0.0\t0.0\t" << _syz << endl;
  os << "0.0\t0.0\t" << _sxz << endl;

  // Invert transformation matrix and ...
  this->Invert();
  // ... and transformation parameters
  this->UpdateParameter();

  return os;
}

irtkCifstream& irtkAffineTransformation::Read(irtkCifstream& from)
{
  int i;
  double data;
  unsigned int magic_no, trans_type, dofs;

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkAffineTransformation::Read: Not a valid transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if ((trans_type != IRTKTRANSFORMATION_AFFINE) && (trans_type != IRTKTRANSFORMATION_RIGID)) {
    cerr << "irtkAffineTransformation::Read: Not a vaild transformation" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&dofs, 1);

  // Initialize rotations and translations
  _tx = _ty = _tz = 0;
  _rx = _ry = _rz = 0;

  // Initialize scale and shears
  _sx  = _sy  = _sz  = 100;
  _sxy = _syz = _sxz = 0;

  // Read data
  for (i = 0; i < int(dofs); i++) {
    from.ReadAsDouble(&data, 1);
    this->Put(i, data);
  }

  return from;
}

irtkCofstream& irtkAffineTransformation::Write(irtkCofstream& to)
{
  int i;

  // Write magic no. for transformations
  unsigned int magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Check if the transformation is really affine
  if ((fabs(_sx - 100.0) < 0.0001) && (fabs(_sy - 100.0) < 0.0001) && (fabs(_sz - 100.0) < 0.0001) &&
      (fabs(_sxy) < 0.0001) && (fabs(_syz) < 0.0001) && (fabs(_sxz) < 0.0001)) {

    // Write transformation type
    unsigned int trans_type = IRTKTRANSFORMATION_RIGID;
    to.WriteAsUInt(&trans_type, 1);

    // Write transformation type
    unsigned int dofs = 6;
    to.WriteAsUInt(&dofs, 1);

    // Write data
    for (i = 0; i < 6; i++) {
      double data = this->Get(i);
      to.WriteAsDouble(&data, 1);
    }

  } else {
    // Write transformation type
    unsigned int trans_type = IRTKTRANSFORMATION_AFFINE;
    to.WriteAsUInt(&trans_type, 1);

    // Write transformation type
    unsigned int dofs = this->NumberOfDOFs();
    to.WriteAsUInt(&dofs, 1);

    // Write data
    for (i = 0; i < this->NumberOfDOFs(); i++) {
      double data = this->Get(i);
      to.WriteAsDouble(&data, 1);
    }
  }

  return to;
}


