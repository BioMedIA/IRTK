
#include <irtkTransformation.h>

#ifdef HAS_VTK

extern bool interactiveVTK;
extern bool displayVTK;
extern bool firstVTK;

extern void update   (irtkPointSet& pset, int nx, int ny, int nz);
extern void visualize(irtkPointSet& pset, int nx, int ny, int nz);

#endif

#define LUTSIZE (double)(BSHLOOKUPTABLESIZE-1)

//double irtkTemporalHomogeneousTransformation::LookupTable   [BSHLOOKUPTABLESIZE][4];
//
//double irtkTemporalHomogeneousTransformation::LookupTable_I [BSHLOOKUPTABLESIZE][4];

void irtkTemporalRigidTransformation::UpdateMatrix()
{
  for (int i = 0; i < _Nt; i++) {
	this->UpdateMatrix(i);
  }
}

void irtkTemporalRigidTransformation::UpdateMatrix(int t)
{
  double tx, ty, tz, rx, ry, rz;
  double cosrx, cosry, cosrz, sinrx, sinry, sinrz;

  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  tx = this->GetTranslationX(t);
  ty = this->GetTranslationY(t);
  tz = this->GetTranslationZ(t);
  rx = this->GetRotationX(t);
  ry = this->GetRotationY(t);
  rz = this->GetRotationZ(t);

  // Update sines and cosines
  cosrx = cos(rx*(M_PI/180.0));
  cosry = cos(ry*(M_PI/180.0));
  cosrz = cos(rz*(M_PI/180.0));
  sinrx = sin(rx*(M_PI/180.0));
  sinry = sin(ry*(M_PI/180.0));
  sinrz = sin(rz*(M_PI/180.0));

  // Add other transformation parameters to transformation matrix
  _matrix[t].Ident();
  _matrix[t](0,0) = cosry*cosrz;
  _matrix[t](0,1) = cosry*sinrz;
  _matrix[t](0,2) = -sinry;
  _matrix[t](0,3) = tx;
  _matrix[t](1,0) = (sinrx*sinry*cosrz-cosrx*sinrz);
  _matrix[t](1,1) = (sinrx*sinry*sinrz+cosrx*cosrz);
  _matrix[t](1,2) = sinrx*cosry;
  _matrix[t](1,3) = ty;
  _matrix[t](2,0) = (cosrx*sinry*cosrz+sinrx*sinrz);
  _matrix[t](2,1) = (cosrx*sinry*sinrz-sinrx*cosrz);
  _matrix[t](2,2) = cosrx*cosry;
  _matrix[t](2,3) = tz;
  _matrix[t](3,3) = 1.0;
}

double irtkTemporalRigidTransformation::GetRotationX(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _rx[t];
}

double irtkTemporalRigidTransformation::GetRotationY(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _ry[t];
}

double irtkTemporalRigidTransformation::GetRotationZ(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _rz[t];
}

double irtkTemporalRigidTransformation::GetTranslationX(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _tx[t];
}

double irtkTemporalRigidTransformation::GetTranslationY(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _ty[t];
}

double irtkTemporalRigidTransformation::GetTranslationZ(int t) const
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  return _tz[t];
}

/// Construct a matrix based on parameters passed in the array.
irtkMatrix irtkTemporalRigidTransformation::Parameters2Matrix(double *params) const
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

void irtkTemporalRigidTransformation::Matrix2Parameters(irtkMatrix m, double* params) const
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

double irtkTemporalRigidTransformation::Get(int i, int t) const
{
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
    default:
      cerr << "irtkTemporalRigidTransformation::Get: No such dof" << endl;
      return 0;
  }
}

double irtkTemporalRigidTransformation::Get(int index, double time)
{
  double t, param;
  double *params = new double [6];

  t = this->TimeToLattice(time);
  if (t < 0 || t >= _Nt)
	cerr << "wrong time "<<t<<endl;

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

    cerr << "irtkTemporalRigidTransformation::Get: No such dof" << endl;
    return 0;
  }
}

double irtkTemporalRigidTransformation::Get(int i) const
{
  int index  = i%6;
  int t = i/6;

  if (t >= _Nt) {
  	cerr << "irtkTemporalRigidTransformation::Get: No such dof" << endl;
  	exit(1);
  }

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
    default:
      cerr << "irtkTemporalRigidTransformation::Get: No such dof" << endl;
      return 0;
  }
}

void irtkTemporalRigidTransformation::Put(int i, int t, double x)
{
  if (t >= _Nt) {
	cerr << "irtkTemporalRigidTransformation::Put: No such dof" << endl;
	exit(1);
  }

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
    default:
      cerr << "irtkTemporalRigidTransformation::Put: No such dof" << endl;
      exit(1);
  }
  this->UpdateMatrix(t);
}

void irtkTemporalRigidTransformation::Put(int i, double x)
{
  int index  = i%6;
  int t = i/6;

  if (t >= _Nt) {
	cerr << "irtkTemporalRigidTransformation::Put: No such dof" << endl;
	exit(1);
  }

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
    default:
      cerr << "irtkTemporalRigidTransformation::Put: No such dof" << endl;
      exit(1);
  }
  this->UpdateMatrix(t);
}

int irtkTemporalRigidTransformation::CheckHeader(char *name)
{
  int n;
  char buffer[255];

  ifstream from(name);

  if (!from) {
    cerr << "irtkTemporalRigidTransformation::CheckHeader: Can't open file "
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
  if (n != 6) {
    return false;
  }

  return true;
}

void irtkTemporalRigidTransformation::Rotate(double& x, double& y, double& z, double t)
{
  double a, b, c;

  this->GetMatrix(t);

  // Pre-multiply point with transformation matrix
  a = _transformMatrix(0, 0)*x+_transformMatrix(0, 1)*y+_transformMatrix(0, 2)*z;
  b = _transformMatrix(1, 0)*x+_transformMatrix(1, 1)*y+_transformMatrix(1, 2)*z;
  c = _transformMatrix(2, 0)*x+_transformMatrix(2, 1)*y+_transformMatrix(2, 2)*z;

  // Copy result back
  x = a;
  y = b;
  z = c;
}

/// B-splines only in the temporal direction, so no after-differentiation
void irtkTemporalRigidTransformation::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double t)
{
  double rx, ry, rz;
  double cosrx, cosry, cosrz, sinrx, sinry, sinrz;
  double *params = new double[6];

  // Get matrix for time point t
  this->GetMatrix(t);

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
      jac[1]  = (cosrx*sinry*cosrz+sinrx*sinrz)*x + (cosrx*sinry*sinrz-sinrx*cosrz)*y + cosrx*cosry*z;
      jac[2]  = (-sinrx*sinry*cosrz+cosrx*sinrz)*x + (-sinrx*sinry*sinrz-cosrx*cosrz)*y - sinrx*cosry*z;
      break;
    case RY:
      // Ensure that the derivatives are expressed with respect to angles not radians
      x *= (M_PI/180.0);
      y *= (M_PI/180.0);
      z *= (M_PI/180.0);
      jac[0]  = -sinry*cosrz*x - sinry*sinrz*y - cosry*z;
      jac[1]  = (sinrx*cosry*cosrz)*x + (sinrx*cosry*sinrz)*y - sinry*sinrx*z;
      jac[2]  = (cosrx*cosry*cosrz)*x + (cosrx*cosry*sinrz)*y - cosrx*sinry*z;
      break;
    case RZ:
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

void irtkTemporalRigidTransformation::Print()
{
  cout.setf(ios::right);
  cout.setf(ios::fixed);
  cout.precision(4);
  for (int i = 0; i < _Nt; i++) {
	cout<<"t: "<<i<<" -> ";
    if (_status[i*6+TX]  == _Active) cout << "tx = " << setw(7) << _tx[i] << " ";
    if (_status[i*6+TY]  == _Active) cout << "ty = " << setw(7) << _ty[i] << " ";
    if (_status[i*6+TZ]  == _Active) cout << "tz = " << setw(7) << _tz[i] << " ";
    if (_status[i*6+RX]  == _Active) cout << "rx = " << setw(7) << _rx[i] << " ";
    if (_status[i*6+RY]  == _Active) cout << "ry = " << setw(7) << _ry[i] << " ";
    if (_status[i*6+RZ]  == _Active) cout << "rz = " << setw(7) << _rz[i] << endl;
  }
  cout.precision(6);
  cout.unsetf(ios::right);
  cout.unsetf(ios::fixed);
}

bool irtkTemporalRigidTransformation::IsIdentity()
{
  bool id = false;
  for (int i = 0; i < _Nt; i++) {
    if ((_tx[i] == 0) && (_ty[i] == 0) && (_tz[i] == 0) &&
    	(_rx[i] == 0) && (_ry[i] == 0) && (_rz[i] == 0)) {
	  id = true;
    } else {
      return false;
    }
  }
  return id;
}

bool irtkTemporalRigidTransformation::IsIdentity(int t)
{
  if ((this->GetTranslationX(t) == 0) && (this->GetTranslationY(t) == 0) && (this->GetTranslationZ(t) == 0) &&
      (this->GetRotationX(t) == 0)    && (this->GetRotationY(t) == 0)    && (this->GetRotationZ(t) == 0)) {
    return true;
  } else {
    return false;
  }
}

void irtkTemporalRigidTransformation::SubdivideOdd()
{
  int i, Nt;
  double * params = new double[6];

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
  tx = this->Allocate(tx, _Nt);
  ty = this->Allocate(ty, _Nt);
  tz = this->Allocate(tz, _Nt);
  rx = this->Allocate(rx, _Nt);
  ry = this->Allocate(ry, _Nt);
  rz = this->Allocate(rz, _Nt);

  // compute new parameters
  for (i = 0; i < _Nt; i++) {
    Matrix2Parameters(_matrix[i], params);
    tx[i] = params[0];
    ty[i] = params[1];
    tz[i] = params[2];
    rx[i] = params[3];
    ry[i] = params[4];
    rz[i] = params[5];
  }

  // Deallocate points
  _tx = this->Deallocate(_tx, Nt);
  _ty = this->Deallocate(_ty, Nt);
  _tz = this->Deallocate(_tz, Nt);
  _rx = this->Deallocate(_rx, Nt);
  _ry = this->Deallocate(_ry, Nt);
  _rz = this->Deallocate(_rz, Nt);
//  delete []_status;

  // Update pointers to control points
  _tx = tx;
  _ty = ty;
  _tz = tz;
  _rx = rx;
  _ry = ry;
  _rz = rz;

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

void irtkTemporalRigidTransformation::SubdivideEven()
{
  int i, Nt;
  double * params = new double[6];

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
  tx = this->Allocate(tx, _Nt);
  ty = this->Allocate(ty, _Nt);
  tz = this->Allocate(tz, _Nt);
  rx = this->Allocate(rx, _Nt);
  ry = this->Allocate(ry, _Nt);
  rz = this->Allocate(rz, _Nt);

  // compute new parameters
  for (i = 0; i < _Nt; i++) {
    Matrix2Parameters(_matrix[i], params);
    tx[i] = params[0];
    ty[i] = params[1];
    tz[i] = params[2];
    rx[i] = params[3];
    ry[i] = params[4];
    rz[i] = params[5];
  }

  // Deallocate points
  _tx = this->Deallocate(_tx, Nt);
  _ty = this->Deallocate(_ty, Nt);
  _tz = this->Deallocate(_tz, Nt);
  _rx = this->Deallocate(_rx, Nt);
  _ry = this->Deallocate(_ry, Nt);
  _rz = this->Deallocate(_rz, Nt);
//  delete []_status;

  // Update pointers to control points
  _tx = tx;
  _ty = ty;
  _tz = tz;
  _rx = rx;
  _ry = ry;
  _rz = rz;

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

istream& irtkTemporalRigidTransformation::Import(istream& is)
{
  int n, t;
  double dummy;
  char buffer[255];

  // Read keyword
  is >> buffer;
  if (strcmp(buffer, "DOF:") != 0) {
    cerr << "irtkTemporalRigidTransformation::Import: Not a valid transformation" << endl;
    exit(1);
  }

  // Read no. of DOFs
  is >> n;
  if (n != 6) {
    cerr << "irtkTemporalRigidTransformation::Import: Not a rigid transformation" << endl;
    exit(1);
  }

  // Read number of time points
  is >> t;
  _Nt = t;

  // Read rigid transformation parameters
  for (int i = 0; i < _Nt; i++) {
    is >> dummy >> dummy >> _tx[i];
    is >> dummy >> dummy >> _ty[i];
    is >> dummy >> dummy >> _tz[i];
    is >> dummy >> dummy >> _rx[i];
    is >> dummy >> dummy >> _ry[i];
    is >> dummy >> dummy >> _rz[i];
  }

  return is;
}

ostream& irtkTemporalRigidTransformation::Export(ostream& os)
{
  // Write rigid transformation parameters
  os << "DOF: 6\n";
  os << _Nt << endl;
  for (int i = 0; i < _Nt; i++) {
    os << "0.0\t0.0\t" << _tx[i] << endl;
    os << "0.0\t0.0\t" << _ty[i] << endl;
    os << "0.0\t0.0\t" << _tz[i] << endl;
    os << "0.0\t0.0\t" << _rx[i] << endl;
    os << "0.0\t0.0\t" << _ry[i] << endl;
    os << "0.0\t0.0\t" << _rz[i] << endl;
  }

  return os;
}

irtkCifstream& irtkTemporalRigidTransformation::Read(irtkCifstream& from)
{
  int i;
  double data, dt, tMin, tMax;
  unsigned int magic_no, trans_type, dofs, Nt;

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkTemporalRigidTransformation::Read: Not a vaild transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if (trans_type != IRTKTRANSFORMATION_RIGID_TEMPORAL) {
    cerr << "irtkTemporalRigidTransformation::Read: Not a vaild rigid transformation " << trans_type << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&dofs, 1);
  if (dofs != 6) {
    cerr << "irtkTemporalRigidTransformation::Read: Invalid no. of dofs = " << dofs << endl;
    exit(1);
  }

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
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    from.ReadAsDouble(&data, 1);
    this->Put(i, data);
  }

  return from;
}

irtkCofstream& irtkTemporalRigidTransformation::Write(irtkCofstream& to)
{
  int i;

  // Write magic no. for transformations
  unsigned int magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

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
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    double data = this->Get(i);
    to.WriteAsDouble(&data, 1);
  }

  return to;
}
