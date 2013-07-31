
#include <irtkTransformation.h>

#define LUTSIZE (double)(BSHLOOKUPTABLESIZE-1)

//double irtkTemporalHomogeneousTransformation::LookupTable   [BSHLOOKUPTABLESIZE][4];
//
//double irtkTemporalHomogeneousTransformation::LookupTable_I [BSHLOOKUPTABLESIZE][4];

void irtkTemporalAffineTransformation::UpdateMatrix()
{
  for (int i = 0; i < _Nt; i++) {
	this->UpdateMatrix(i);
  }
}

void irtkTemporalAffineTransformation::UpdateMatrix(int t)
{
  double sx, sy, sz, sxy, syz, sxz;

  // Update rigid transformation
  irtkTemporalRigidTransformation::UpdateMatrix(t);

  sx = this->GetScaleX(t);
  sy = this->GetScaleY(t);
  sz = this->GetScaleZ(t);

  sxy = this->GetShearXY(t);
  syz = this->GetShearYZ(t);
  sxz = this->GetShearXZ(t);

  // Update affine transformation: Add shearing
  irtkMatrix skew(4, 4);
  skew.Ident();
  skew(0, 1) = tan(sxy*(M_PI/180.0));
  skew(0, 2) = tan(sxz*(M_PI/180.0));
  skew(1, 2) = tan(syz*(M_PI/180.0));
  _matrix[t] *= skew;

  // Update affine transformation: Add scaling
  irtkMatrix scale(4, 4);
  scale.Ident();
  scale(0, 0) = sx / 100.0;
  scale(1, 1) = sy / 100.0;
  scale(2, 2) = sz / 100.0;
  _matrix[t] *= scale;
}

double irtkTemporalAffineTransformation::GetScaleX(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _sx[t];
}

double irtkTemporalAffineTransformation::GetScaleY(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _sy[t];
}

double irtkTemporalAffineTransformation::GetScaleZ(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _sz[t];
}

double irtkTemporalAffineTransformation::GetShearXY(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _sxy[t];
}

double irtkTemporalAffineTransformation::GetShearYZ(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _syz[t];
}

double irtkTemporalAffineTransformation::GetShearXZ(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _sxz[t];
}

/// Construct a matrix based on parameters passed in the array.
irtkMatrix irtkTemporalAffineTransformation::Parameters2Matrix(double *params) const
{
  // Get the rigid components.
  irtkMatrix matrix;
  matrix = this->irtkTemporalRigidTransformation::Parameters2Matrix(params);

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

//void irtkTemporalAffineTransformation::Matrix2Parameters(irtkMatrix m, double* params) const
//{
//  int i, j;
//  irtkVector Z(3);
//  irtkMatrix U(3, 3);
//  irtkMatrix V(3, 3);
//  irtkMatrix Tmp(3, 3);
//
//  // compute scaling
//  for (i = 0; i < 3; i++) {
//	for (j = 0; j < 3; j++) {
//	  Tmp(i, j) = m(i, j);
//	}
//  }
//
//  // SVD = Rotation*Scaling*Rotation
//  Tmp.SVD(U, Z, V);
//
//  params[SX] = Z(0);
//  params[SY] = Z(1);
//  params[SZ] = Z(2);
//
//  // remove scaling from matrix
//  for (i = 0; i < 3; i++) {
//	for (j = 0; j < 4; j++) {
//	  m(i, j) = m(i, j)/params[6+i];
//	}
//  }
//
//  this->irtkTemporalRigidTransformation::Matrix2Parameters(m, params);
//}

// Use current matrix to extract the transformation parameters based on 12
// DOF model, see:
//    http://www.acm.org/pubs/tog/GraphicsGems/gemsii/unmatrix.c
// (from Graphics Gems II)
void irtkTemporalAffineTransformation::Matrix2Parameters(irtkMatrix matrix, double* params) const
{
  // Use _matrix to evaluate the affine transformation parameters.
  // It is assumed that there is no perspective transformation.
  int i;
  double TOL = 0.000001;
  double tansxy, tansxz, tansyz;

  if (fabs(matrix(3,3) - 1.0) > TOL) {
    cerr << "irtkAffineTransformation::UpdateParameter" << endl;
    cerr << "Value at matrix(3,3) must equal 1." << endl;
    exit(1);
  }

  if (fabs(matrix.Det()) < TOL) {
    cerr << "irtkAffineTransformation::UpdateParameter" << endl;
    cerr << "Matrix singular (or very close to singular!)." << endl;
    exit(1);
  }

  // First Part Of Graphics Gems Code Ignored Because It Relates To
  // Perspective Transformation.
  if (fabs(matrix(3, 0)) > TOL ||
      fabs(matrix(3, 1)) > TOL ||
      fabs(matrix(3, 2)) > TOL ) {
    cerr << "irtkAffineTransformation::UpdateParameter" << endl;
    cerr << "Matrix contains perspective components." << endl;
    exit(1);
  }

  irtkMatrix copy(4, 4);
  copy = matrix;

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
  params[SX] = col_0.Norm();
  col_0 /= params[SX];

  // Compute XY shear factor and make 2nd col orthogonal to 1st.
  tansxy = col_0.ScalarProduct(col_1);
  col_1 = col_1 - col_0 * tansxy;

  // Actually, tansxy and col_1 are still to large by a factor of sy.
  // Now, compute Y scale and normalize 2nd col and rescale tansxy.
  params[SY] = col_1.Norm();
  col_1  /= params[SY];
  tansxy /= params[SY];

  // Compute XZ and YZ shears, orthogonalize 3rd col
  tansxz = col_0.ScalarProduct(col_2);
  col_2 = col_2 - col_0 * tansxz;

  tansyz = col_1.ScalarProduct(col_2);
  col_2 = col_2 - col_1 * tansyz;

  // Actually, tansxz, tansyz and col_2 are still to large by a factor of
  // sz.  Next, get Z scale, normalize 3rd col and scale tansxz and tansyz.
  params[SZ] = col_2.Norm();
  col_2  /= params[SZ];
  tansxz /= params[SZ];
  tansyz /= params[SZ];

  // At this point, the columns are orthonormal.  Check for a coordinate
  // system flip.  If the determinant is -1, then negate the matrix and the
  // scaling factors.
  irtkVector col_1_x_col_2;
  col_1_x_col_2.Initialize(3);
  col_1_x_col_2 = col_1.CrossProduct(col_2);

  if (col_0.ScalarProduct(col_1_x_col_2) < 0) {
	params[SX] *= -1;
	params[SY] *= -1;
	params[SZ] *= -1;
    col_0 *= -1;
    col_1 *= -1;
    col_2 *= -1;
  }

  // Retrieve the shear angles in degrees.
  params[SXY] = atan(tansxy) * 180.0 / M_PI ;
  params[SXZ] = atan(tansxz) * 180.0 / M_PI ;
  params[SYZ] = atan(tansyz) * 180.0 / M_PI ;

  // Convert scales to percentages.
  params[SX] *= 100;
  params[SY] *= 100;
  params[SZ] *= 100;

  // Now get the rigid transformation parameters.
  double rigidParams[6];
  // Put the rotation matrix components into the upper left 3x3 submatrix.
  for (i = 0; i < 3; ++i) {
    copy(i, 0) = col_0(i);
    copy(i, 1) = col_1(i);
    copy(i, 2) = col_2(i);
  }

  this->irtkTemporalRigidTransformation::Matrix2Parameters(copy, rigidParams);

  params[TX] = rigidParams[TX];
  params[TY] = rigidParams[TY];
  params[TZ] = rigidParams[TZ];
  params[RX] = rigidParams[RX];
  params[RY] = rigidParams[RY];
  params[RZ] = rigidParams[RZ];
}

double irtkTemporalAffineTransformation::Get(int i, int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  switch (i) {
    case TX:
      return _tx[t];
      break;
    case TY:
      return _ty[t];
      break;
    case TZ:
      return _tz[t];
      break;
    case RX:
      return _rx[t];
      break;
    case RY:
      return _ry[t];
      break;
    case RZ:
      return _rz[t];
      break;
    case SX:
      return _sx[t];
      break;
    case SY:
      return _sy[t];
      break;
    case SZ:
      return _sz[t];
      break;
    case SXY:
      return _sxy[t];
      break;
    case SYZ:
      return _syz[t];
      break;
    case SXZ:
      return _sxz[t];
      break;
    default:
      cerr << "irtkTemporalAffineTransformation::Get: No such dof" << endl;
      return 0;
  }
}

double irtkTemporalAffineTransformation::Get(int index, double time)
{
  double t, param;
  double *params = new double [12];

  t = this->TimeToLattice(time);
  if (t < 0 || t >= _Nt)
	cerr << "wrong time "<<t<<" (out of bounds)"<<endl;

  this->GetMatrix(time);
  this->Matrix2Parameters(_transformMatrix, params);

  if (index < this->NumberOfParams()) {
	param = params[index];

	delete params;
	params = NULL;

	return param;
  } else {
	delete params;
	params = NULL;

    cerr << "irtkTemporalAffineTransformation::Get: No such dof" << endl;
    return 0;
  }
}

double irtkTemporalAffineTransformation::Get(int i) const
{
  int index  = i%12;
  int t = i/12;

  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  switch (index) {
    case 0:
      return _tx[t];
      break;
    case 1:
      return _ty[t];
      break;
    case 2:
      return _tz[t];
      break;
    case 3:
      return _rx[t];
      break;
    case 4:
      return _ry[t];
      break;
    case 5:
      return _rz[t];
      break;
    case 6:
      return _sx[t];
      break;
    case 7:
      return _sy[t];
      break;
    case 8:
      return _sz[t];
      break;
    case 9:
      return _sxy[t];
      break;
    case 10:
      return _syz[t];
      break;
    case 11:
      return _sxz[t];
      break;
    default:
      cerr << "irtkTemporalAffineTransformation::Get: No such dof" << endl;
      return 0;
  }
}

void irtkTemporalAffineTransformation::Put(int i, int t, double x)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  switch (i) {
    case TX:
      _tx[t] = x;
      break;
    case TY:
      _ty[t] = x;
      break;
    case TZ:
      _tz[t] = x;
      break;
    case RX:
      _rx[t] = x;
      break;
    case RY:
      _ry[t] = x;
      break;
    case RZ:
      _rz[t] = x;
      break;
    case SX:
      _sx[t] = x;
      break;
    case SY:
      _sy[t] = x;
      break;
    case SZ:
      _sz[t] = x;
      break;
    case SXY:
      _sxy[t] = x;
      break;
    case SYZ:
      _syz[t] = x;
      break;
    case SXZ:
      _sxz[t] = x;
      break;
    default:
      cerr << "irtkTemporalAffineTransformation::Put(): No such dof" << endl;
      exit(1);
  }
  this->UpdateMatrix(t);
}

void irtkTemporalAffineTransformation::Put(int i, double x)
{
  int index  = i%12;
  int t = i/12;

  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  switch (index) {
    case 0:
      _tx[t] = x;
      break;
    case 1:
      _ty[t] = x;
      break;
    case 2:
      _tz[t] = x;
      break;
    case 3:
      _rx[t] = x;
      break;
    case 4:
      _ry[t] = x;
      break;
    case 5:
      _rz[t] = x;
      break;
    case 6:
      _sx[t] = x;
      break;
    case 7:
      _sy[t] = x;
      break;
    case 8:
      _sz[t] = x;
      break;
    case 9:
      _sxy[t] = x;
      break;
    case 10:
      _syz[t] = x;
      break;
    case 11:
      _sxz[t] = x;
      break;
    default:
      cerr << "irtkTemporalAffineTransformation::Put(): No such dof" << endl;
      exit(1);
  }
  this->UpdateMatrix(t);
}

/// B-splines only in the temporal direction, so no after-differentiation
void irtkTemporalAffineTransformation::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double t)
{
  double tmp;
  double rx, ry, rz;
  double cosrx, cosry, cosrz, sinrx, sinry, sinrz;
  double sx, sy, sz, sxy, sxz, syz;
  double tansxy, tansxz, tansyz;
  double *params = new double[12];

  // Get matrix for time point t
  this->GetMatrix(t);
  _tCurrent = t;

  // compute parameters for this transformation matrix
  this->Matrix2Parameters(_transformMatrix, params);

  rx = params[RX];
  ry = params[RY];
  rz = params[RZ];
  cosrx = cos(rx*(M_PI/180.0));
  cosry = cos(ry*(M_PI/180.0));
  cosrz = cos(rz*(M_PI/180.0));
  sinrx = sin(rx*(M_PI/180.0));
  sinry = sin(ry*(M_PI/180.0));
  sinrz = sin(rz*(M_PI/180.0));

  // Rotation matrix
  irtkMatrix r(3, 3);
  r(0, 0) = cosry*cosrz;
  r(0, 1) = cosry*sinrz;
  r(0, 2) = -sinry;
  r(1, 0) = (sinrx*sinry*cosrz-cosrx*sinrz);
  r(1, 1) = (sinrx*sinry*sinrz+cosrx*cosrz);
  r(1, 2) = sinrx*cosry;
  r(2, 0) = (cosrx*sinry*cosrz+sinrx*sinrz);
  r(2, 1) = (cosrx*sinry*sinrz-sinrx*cosrz);
  r(2, 2) = cosrx*cosry;

  sx = params[SX];
  sy = params[SY];
  sz = params[SZ];
  sxy = params[SXY];
  sxz = params[SXZ];
  syz = params[SYZ];
  tansxy = tan(sxy*(M_PI/180.0));
  tansxz = tan(sxz*(M_PI/180.0));
  tansyz = tan(syz*(M_PI/180.0));

  switch (dof) {
    case SX:
      jac[0] = r(0, 0) * x / 100.0;
      jac[1] = r(1, 0) * x / 100.0;
      jac[2] = r(2, 0) * x / 100.0;
      break;
    case SY:
      jac[0] = r(0, 0) * tansxy * y / 100.0 + r(0, 1) * y / 100.0;
      jac[1] = r(1, 0) * tansxy * y / 100.0 + r(1, 1) * y / 100.0;
      jac[2] = r(2, 0) * tansxy * y / 100.0 + r(2, 1) * y / 100.0;
      break;
    case SZ:
      jac[0] = r(0, 0) * tansxz * z / 100.0 + r(0, 1) * tansyz * z / 100.0 + r(0, 2) * z / 100.0;
      jac[1] = r(1, 0) * tansxz * z / 100.0 + r(1, 1) * tansyz * z / 100.0 + r(1, 2) * z / 100.0;
      jac[2] = r(2, 0) * tansxz * z / 100.0 + r(2, 1) * tansyz * z / 100.0 + r(2, 2) * z / 100.0;
      break;
    case SXY:
      tmp = (1.0 / pow(cos(sxy*(M_PI/180.0)), 2.0))*(M_PI/180.0) * sy * y / 100.0;
      jac[0] = r(0, 0) * tmp;
      jac[1] = r(1, 0) * tmp;
      jac[2] = r(2, 0) * tmp;
      break;
    case SYZ:
      tmp = (1.0 / pow(cos(syz*(M_PI/180.0)), 2.0))*(M_PI/180.0) * sz * z / 100.0;
      jac[0] = r(0, 1) * tmp;
      jac[1] = r(1, 1) * tmp;
      jac[2] = r(2, 1) * tmp;
      break;
    case SXZ:
      tmp = (1.0 / pow(cos(sxz*(M_PI/180.0)), 2.0))*(M_PI/180.0) * sz * z / 100.0;
      jac[0] = r(0, 0) * tmp;
      jac[1] = r(1, 0) * tmp;
      jac[2] = r(2, 0) * tmp;
      break;
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
      x = sx / 100.0 * x + tansxy * y * sy / 100.0 + tansxz * z * sz / 100.0;
      y = sy / 100.0 * y + tansyz * z * sz / 100.0;
      z = sz / 100.0 * z;
      // Ensure that the derivatives are expressed with respect to angles not radians
      x *= (M_PI/180.0);
      y *= (M_PI/180.0);
      z *= (M_PI/180.0);
      jac[0]  = 0;
      jac[1]  = (cosrx*sinry*cosrz+sinrx*sinrz)*x + (cosrx*sinry*sinrz-sinrx*cosrz)*y + cosrx*cosry*z;
      jac[2]  = (-sinrx*sinry*cosrz+cosrx*sinrz)*x + (-sinrx*sinry*sinrz-cosrx*cosrz)*y - sinrx*cosry*z;
      break;
    case RY:
      x = sx / 100.0 * x + tansxy * y * sy / 100.0 + tansxz * z * sz / 100.0;
      y = sy / 100.0 * y + tansyz * z * sz / 100.0;
      z = sz / 100.0 * z;
      // Ensure that the derivatives are expressed with respect to angles not radians
      x *= (M_PI/180.0);
      y *= (M_PI/180.0);
      z *= (M_PI/180.0);
      jac[0]  = -sinry*cosrz*x - sinry*sinrz*y - cosry*z;
      jac[1]  = (sinrx*cosry*cosrz)*x + (sinrx*cosry*sinrz)*y - sinry*sinrx*z;
      jac[2]  = (cosrx*cosry*cosrz)*x + (cosrx*cosry*sinrz)*y - cosrx*sinry*z;
      break;
    case RZ:
      x = sx / 100.0 * x + tansxy * y * sy / 100.0 + tansxz * z * sz / 100.0;
      y = sy / 100.0 * y + tansyz * z * sz / 100.0;
      z = sz / 100.0 * z;
      // Ensure that the derivatives are expressed with respect to angles not radians
      x *= (M_PI/180.0);
      y *= (M_PI/180.0);
      z *= (M_PI/180.0);
      jac[0]  = -sinrz*cosry*x + cosry*cosrz*y;
      jac[1]  = (-sinrx*sinry*sinrz-cosrx*cosrz)*x + (sinrx*sinry*cosrz-cosrx*sinrz)*y;
      jac[2]  = (-cosrx*sinry*sinrz+sinrx*cosrz)*x + (cosrx*sinry*cosrz+sinrx*sinrz)*y;
      break;
    default:
      cerr << this->NameOfClass() << "::JacobianDOFs(): No such dof = " << dof << endl;
      exit(1);
      break;
  }

  delete params;
  params = NULL;
}

int irtkTemporalAffineTransformation::CheckHeader(char *name)
{
  int n;
  char buffer[255];

  ifstream from(name);

  if (!from) {
    cerr << "irtkTemporalAffineTransformation::CheckHeader: Can't open file "
    << name << "\n";
    exit(1);
  }

  // Read keyword
  from >> buffer;
  if (strcmp(buffer, "DOF:") != 0) {
    return false;
  }

  // Read no. of DOFs
  from >> n;
  if ((n != 6) && (n != 9) && (n != 12) && (n != 15)) {
    return false;
  }

  return true;
}

void irtkTemporalAffineTransformation::Print()
{
  int nDOF = this->NumberOfDOFs()/_Nt;

  cout.setf(ios::right);
  cout.setf(ios::fixed);
  cout.precision(4);
  for (int i = 0; i < _Nt; i++) {
	cout<<"t: "<<i<<" -> ";
    if (_status[i*nDOF+TX]  == _Active) cout << "tx = "  << setw(7) << _tx[i] << " ";
    if (_status[i*nDOF+TY]  == _Active) cout << "ty = "  << setw(7) << _ty[i] << " ";
    if (_status[i*nDOF+TZ]  == _Active) cout << "tz = "  << setw(7) << _tz[i] << " ";
    if (_status[i*nDOF+RX]  == _Active) cout << "rx = "  << setw(7) << _rx[i] << " ";
    if (_status[i*nDOF+RY]  == _Active) cout << "ry = "  << setw(7) << _ry[i] << " ";
    if (_status[i*nDOF+RZ]  == _Active) cout << "rz = "  << setw(7) << _rz[i] << " ";
//    cout << endl;
    if (_status[i*nDOF+SX]  == _Active) cout << "sx = "  << setw(7) << _sx[i] << " ";
    if (_status[i*nDOF+SY]  == _Active) cout << "sy = "  << setw(7) << _sy[i] << " ";
    if (_status[i*nDOF+SZ]  == _Active) cout << "sz = "  << setw(7) << _sz[i] << " ";
    if (_status[i*nDOF+SXY] == _Active) cout << "sxy = " << setw(7) << _sxy[i] << " ";
    if (_status[i*nDOF+SYZ] == _Active) cout << "syz = " << setw(7) << _syz[i] << " ";
    if (_status[i*nDOF+SXZ] == _Active) cout << "sxz = " << setw(7) << _sxz[i] << " ";
    cout << endl;
  }

  cout.precision(6);
  cout.unsetf(ios::right);
  cout.unsetf(ios::fixed);
}

bool irtkTemporalAffineTransformation::IsIdentity()
{
  bool id = false;
  for (int i = 0; i < _Nt; i++) {
	if ((_tx[i] == 0)  && (_ty[i] == 0)  && (_tz[i] == 0)  &&
		(_rx[i] == 0)  && (_ry[i] == 0)  && (_rz[i] == 0)  &&
		(_sx[i] == 100.0)  && (_sy[i] == 100.0)  && (_sz[i] == 100.0)  &&
		(_sxy[i] == 0) && (_syz[i] == 0) && (_sxz[i] == 0)) {
	  id = true;
	} else {
	  return false;
	}
  }
  return id;
}

bool irtkTemporalAffineTransformation::IsIdentity(int t)
{
  if ((this->GetTranslationX(t) == 0) && (this->GetTranslationY(t) == 0) && (this->GetTranslationZ(t) == 0) &&
	  (this->GetRotationX(t) == 0)    && (this->GetRotationY(t) == 0)    && (this->GetRotationZ(t) == 0)    &&
      (this->GetScaleX(t) == 100.0)   && (this->GetScaleY(t) == 100.0)   && (this->GetScaleZ(t) == 100.0)   &&
      (this->GetShearXY(t) == 0)      && (this->GetShearXZ(t) == 0)      && (this->GetShearYZ(t) == 0)) {
    return true;
  } else {
    return false;
  }
}

void irtkTemporalAffineTransformation::SubdivideOdd()
{
  int i, Nt;
  double * params = new double[12];

  // save Nt because it will be overwriten
  Nt = _Nt;

  // Subdivide matrizes
  irtkTemporalHomogeneousTransformation::SubdivideOdd();

  // Allocate memory for new control points
  double * tx = NULL;
  double * ty = NULL;
  double * tz = NULL;
  double * rx = NULL;
  double * ry = NULL;
  double * rz = NULL;
  double * sx = NULL;
  double * sy = NULL;
  double * sz = NULL;
  double * sxy = NULL;
  double * sxz = NULL;
  double * syz = NULL;
  tx = this->Allocate(tx, _Nt);
  ty = this->Allocate(ty, _Nt);
  tz = this->Allocate(tz, _Nt);
  rx = this->Allocate(rx, _Nt);
  ry = this->Allocate(ry, _Nt);
  rz = this->Allocate(rz, _Nt);
  sx = this->Allocate(sx, _Nt);
  sy = this->Allocate(sy, _Nt);
  sz = this->Allocate(sz, _Nt);
  sxy = this->Allocate(sxy, _Nt);
  sxz = this->Allocate(sxz, _Nt);
  syz = this->Allocate(syz, _Nt);

  // compute new parameters
  for (i = 0; i < _Nt; i++) {
    Matrix2Parameters(_matrix[i], params);
    tx[i] = params[0];
    ty[i] = params[1];
    tz[i] = params[2];
    rx[i] = params[3];
    ry[i] = params[4];
    rz[i] = params[5];
    sx[i] = params[6];
    sy[i] = params[7];
    sz[i] = params[8];
    sxy[i] = params[9];
    syz[i] = params[10];
    sxz[i] = params[11];
  }

  // Deallocate points
  _tx = this->Deallocate(_tx, Nt);
  _ty = this->Deallocate(_ty, Nt);
  _tz = this->Deallocate(_tz, Nt);
  _rx = this->Deallocate(_rx, Nt);
  _ry = this->Deallocate(_ry, Nt);
  _rz = this->Deallocate(_rz, Nt);
  _sx = this->Deallocate(_sx, Nt);
  _sy = this->Deallocate(_sy, Nt);
  _sz = this->Deallocate(_sz, Nt);
  _sxy = this->Deallocate(_sxy, Nt);
  _sxz = this->Deallocate(_sxz, Nt);
  _syz = this->Deallocate(_syz, Nt);
//  delete []_status;

  // Update pointers to control points
  _tx = tx;
  _ty = ty;
  _tz = tz;
  _rx = rx;
  _ry = ry;
  _rz = rz;
  _sx = sx;
  _sy = sy;
  _sz = sz;
  _sxy = sxy;
  _sxz = sxz;
  _syz = syz;

  // Increase number of control points --> done in Homogeneous Subdivision
//  _Nt = 2*_Nt-1;

  // Recalculate control point spacing --> done in Homogeneous Subdivision
//  _dt = _dt/2;

  // Initialize memory for control point status --> done in Homogeneous Subdivision
//  _status = new _Status[this->NumberOfDOFs()];
//  for (i = 0; i < this->NumberOfDOFs(); i++) {
//    _status[i] = _Active;
//  }

  delete params;
  params = NULL;
}

void irtkTemporalAffineTransformation::SubdivideEven()
{
  int i, Nt;
  double * params = new double[12];

  // save Nt because it will be overwriten
  Nt = _Nt;

  // Subdivide matrizes
  irtkTemporalHomogeneousTransformation::SubdivideEven();

  // Allocate memory for new control points
  double * tx = NULL;
  double * ty = NULL;
  double * tz = NULL;
  double * rx = NULL;
  double * ry = NULL;
  double * rz = NULL;
  double * sx = NULL;
  double * sy = NULL;
  double * sz = NULL;
  double * sxy = NULL;
  double * sxz = NULL;
  double * syz = NULL;
  tx = this->Allocate(tx, _Nt);
  ty = this->Allocate(ty, _Nt);
  tz = this->Allocate(tz, _Nt);
  rx = this->Allocate(rx, _Nt);
  ry = this->Allocate(ry, _Nt);
  rz = this->Allocate(rz, _Nt);
  sx = this->Allocate(sx, _Nt);
  sy = this->Allocate(sy, _Nt);
  sz = this->Allocate(sz, _Nt);
  sxy = this->Allocate(sxy, _Nt);
  sxz = this->Allocate(sxz, _Nt);
  syz = this->Allocate(syz, _Nt);

  // compute new parameters
  for (i = 0; i < _Nt; i++) {
    Matrix2Parameters(_matrix[i], params);
    tx[i] = params[0];
    ty[i] = params[1];
    tz[i] = params[2];
    rx[i] = params[3];
    ry[i] = params[4];
    rz[i] = params[5];
    sx[i] = params[6];
    sy[i] = params[7];
    sz[i] = params[8];
    sxy[i] = params[9];
    syz[i] = params[10];
    sxz[i] = params[11];
  }

  // Deallocate points
  _tx = this->Deallocate(_tx, Nt);
  _ty = this->Deallocate(_ty, Nt);
  _tz = this->Deallocate(_tz, Nt);
  _rx = this->Deallocate(_rx, Nt);
  _ry = this->Deallocate(_ry, Nt);
  _rz = this->Deallocate(_rz, Nt);
  _sx = this->Deallocate(_sx, Nt);
  _sy = this->Deallocate(_sy, Nt);
  _sz = this->Deallocate(_sz, Nt);
  _sxy = this->Deallocate(_sxy, Nt);
  _sxz = this->Deallocate(_sxz, Nt);
  _syz = this->Deallocate(_syz, Nt);
//  delete []_status;

  // Update pointers to control points
  _tx = tx;
  _ty = ty;
  _tz = tz;
  _rx = rx;
  _ry = ry;
  _rz = rz;
  _sx = sx;
  _sy = sy;
  _sz = sz;
  _sxy = sxy;
  _sxz = sxz;
  _syz = syz;

  // Increase number of control points --> done in Homogeneous Subdivision
//  _Nt = 2*_Nt-1;

  // Recalculate control point spacing --> done in Homogeneous Subdivision
//  _dt = _dt/2;

  // Initialize memory for control point status --> done in Homogeneous Subdivision
//  _status = new _Status[this->NumberOfDOFs()];
//  for (i = 0; i < this->NumberOfDOFs(); i++) {
//    _status[i] = _Active;
//  }

  delete params;
  params = NULL;
}

istream& irtkTemporalAffineTransformation::Import(istream& is)
{
  int n, t;
  double dummy;
  char buffer[255];
  double params[15];
  memset ((void*)&params[0], 0, 15 * sizeof(double));

  // Read keyword
  is >> buffer;
  if (strcmp(buffer, "DOF:") != 0) {
    cerr << "irtkTemporalAffineTransformation::operator>>: Not a valid transformation"
    << endl;
    exit(1);
  }

  // Read no. of DOFs
  is >> n;
  if ((n != 6)  && (n != 9) && (n != 12) && (n != 15)) {
    cerr << "irtkTemporalAffineTransformation::operator>>: Not an affine or rigid "
    << "transformation" << endl;
    exit(1);
  }

  // Read number of time points
  is >> t;
  _Nt = t;

  // Free memory for control points if necessary
  if (_tx != NULL) _tx = this->Deallocate(_tx, _Nt);
  if (_ty != NULL) _ty = this->Deallocate(_ty, _Nt);
  if (_tz != NULL) _tz = this->Deallocate(_tz, _Nt);
  if (_rx != NULL) _rx = this->Deallocate(_rx, _Nt);
  if (_ry != NULL) _ry = this->Deallocate(_ry, _Nt);
  if (_rz != NULL) _rz = this->Deallocate(_rz, _Nt);
  if (_sx != NULL) _sx = this->Deallocate(_sx, _Nt);
  if (_sy != NULL) _sy = this->Deallocate(_sy, _Nt);
  if (_sz != NULL) _sz = this->Deallocate(_sz, _Nt);
  if (_sxy != NULL) _sxy = this->Deallocate(_sxy, _Nt);
  if (_sxz != NULL) _sxz = this->Deallocate(_sxz, _Nt);
  if (_syz != NULL) _syz = this->Deallocate(_syz, _Nt);

  // Allocate parameters
  _tx = this->Allocate(_tx, _Nt);
  _ty = this->Allocate(_ty, _Nt);
  _tz = this->Allocate(_tz, _Nt);
  _rx = this->Allocate(_rx, _Nt);
  _ry = this->Allocate(_ry, _Nt);
  _rz = this->Allocate(_rz, _Nt);
  _sx = this->Allocate(_sx, _Nt);
  _sy = this->Allocate(_sy, _Nt);
  _sz = this->Allocate(_sz, _Nt);
  _sxy = this->Allocate(_sxy, _Nt);
  _sxz = this->Allocate(_sxz, _Nt);
  _syz = this->Allocate(_syz, _Nt);

  // Read rigid transformation parameters
  for (int i = 0; i < _Nt; i++) {
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
      this->SetParameters(&params[0], i);
    } else { // n == 15
      cerr << "irtkTemporalAffineTransformation::Import: 15DOF not implemented" << endl;
      exit(1);
    }
  }

  return is;
}

ostream& irtkTemporalAffineTransformation::Export(ostream& os)
{
  // Write affine transformation parameters
  os << "DOF: " << this->NumberOfDOFs()/_Nt << endl;
  os << _Nt << endl;
  for (int i = 0; i < _Nt; i++) {
    os << "0.0\t0.0\t" << _tx[i] << endl;
    os << "0.0\t0.0\t" << _ty[i] << endl;
    os << "0.0\t0.0\t" << _tz[i] << endl;
    os << "0.0\t0.0\t" << _rx[i] << endl;
    os << "0.0\t0.0\t" << _ry[i] << endl;
    os << "0.0\t0.0\t" << _rz[i] << endl;
    os << "0.0\t0.0\t" << _sx[i]/100.0 << endl;
    os << "0.0\t0.0\t" << _sy[i]/100.0 << endl;
    os << "0.0\t0.0\t" << _sz[i]/100.0 << endl;
    os << "0.0\t0.0\t" << _sxy[i] << endl;
    os << "0.0\t0.0\t" << _syz[i] << endl;
    os << "0.0\t0.0\t" << _sxz[i] << endl;
  }

  return os;
}

irtkCifstream& irtkTemporalAffineTransformation::Read(irtkCifstream& from)
{
  int i;
  double data, dt, tMin, tMax;
  unsigned int magic_no, trans_type, dofs, Nt;

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkTemporalAffineTransformation::Read: Not a valid transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if ((trans_type != IRTKTRANSFORMATION_AFFINE_TEMPORAL) && (trans_type != IRTKTRANSFORMATION_RIGID_TEMPORAL)) {
    cerr << "irtkTemporalAffineTransformation::Read: Not a vaild transformation" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&dofs, 1);

  // Read number of time points
  from.ReadAsUInt(&Nt, 1);
  _Nt = Nt;

  // Read control point spacing
  from.ReadAsDouble(&dt, 1);
  _dt = dt;

  // Read interval of time values for which the transformation is defined
  from.ReadAsDouble(&tMin, 1);
  _tMin = tMin;
  from.ReadAsDouble(&tMax, 1);
  _tMax = tMax;

  this->Initialize(dt, tMin, tMax);

  // Read data
  for (i = 0; i < int(dofs*_Nt); i++) {
    from.ReadAsDouble(&data, 1);
    this->Put(i, data);
  }

  return from;
}

irtkCofstream& irtkTemporalAffineTransformation::Write(irtkCofstream& to)
{
  int i, t, a_r = 1;

  // Write magic no. for transformations
  unsigned int magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Check if the transformation is really affine
  for (t = 0; t < _Nt; t++) {
    if ((fabs(_sx[t] - 100.0) >= 0.0001) && (fabs(_sy[t] - 100.0) >= 0.0001) && (fabs(_sz[t] - 100.0) >= 0.0001) &&
        (fabs(_sxy[t]) >= 0.0001) && (fabs(_syz[t]) >= 0.0001) && (fabs(_sxz[t]) >= 0.0001)) {
    	a_r = 0;
    }
  }
  a_r=0;
  if (a_r==1) {
    // Write transformation type
    unsigned int trans_type = IRTKTRANSFORMATION_RIGID_TEMPORAL;
    to.WriteAsUInt(&trans_type, 1);

    // Write transformation type
    unsigned int dofs = 6;
    to.WriteAsUInt(&dofs, 1);

    // Write number of time points
    unsigned int Nt = _Nt;
    to.WriteAsUInt(&Nt, 1);

    // Write control point spacing
    double dt = _dt;
    to.WriteAsDouble(&dt, 1);

    // Write interval of time values for which the transformation is defined
    double tMin = _tMin;
    to.WriteAsDouble(&tMin, 1);
    double tMax = _tMax;
    to.WriteAsDouble(&tMax, 1);

    // Write data
    for (i = 0; i < 6*this->_Nt; i++) {
      double data = this->Get(i);
      to.WriteAsDouble(&data, 1);
    }
  } else {
    // Write transformation type
    unsigned int trans_type = IRTKTRANSFORMATION_AFFINE_TEMPORAL;
    to.WriteAsUInt(&trans_type, 1);

    // Write transformation type
    unsigned int dofs = this->NumberOfDOFs()/_Nt;
    to.WriteAsUInt(&dofs, 1);

    // Write number of time points
    unsigned int Nt = this->_Nt;
    to.WriteAsUInt(&Nt, 1);

    // Write control point spacing
    double dt = _dt;
    to.WriteAsDouble(&dt, 1);

    // Write interval of time values for which the transformation is defined
    double tMin = _tMin;
    to.WriteAsDouble(&tMin, 1);
    double tMax = _tMax;
    to.WriteAsDouble(&tMax, 1);

    // Write data
    for (i = 0; i < this->NumberOfDOFs(); i++) {
      double data = this->Get(i);
      to.WriteAsDouble(&data, 1);
    }
  }

  return to;
}


