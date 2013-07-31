
#include <irtkTransformation.h>

inline double cspline(double x)
{
  double xi, xii, xiii;

  xi  = fabs(x);
  xii = xi*xi;
  xiii = xii*xi;
  if (xi < 1) {
    return (xiii - 2*xii + 1);
  } else {
    if (xi < 2) {
      return (-xiii + 5*xii - 8*xi + 4);
    } else {
      return 0;
    }
  }
}

irtkMatrix irtkTemporalHomogeneousTransformation::GetMatrix(double time)
{
  double t;

  t = this->TimeToLattice(time);
//cout<<"irtkTemporalHomogeneousTransformation::GetMatrix t = "<<t<<" (time = "<<time<<")"<<endl;
//  // Check if matrix for this time is already computed
//  if (_tCurrent == time) {
//	return _transformMatrix;
//  }

  // Check if there is some work to do
  if ((t < -2) || (t > _Nt+1)) {
	_tCurrent = time;
	_transformMatrix.Ident();
	return _transformMatrix;
  }

  if (_mode == CSpline) {

	int l, d;
	double *weight = new double[4];
	irtkMatrix *M = new irtkMatrix[4];

	// compute weights for interpolation
    d = (int)floor(t)-1;

	for (l = 0; l < 4; l++) {
	  weight[l] = cspline(l+d-t);
	}
	// get matrizes for interpolation
	for (l = 0; l < 4; l++) {
	  M[l] = _matrix[d+l];
	}

	// compute new matrix
	_tCurrent = time;
	_transformMatrix = FrechetMean(M, weight, 4);

	delete[] M;
	M = NULL;
	delete weight;
	weight = NULL;
  } else {
	if (_mode == Linear) {
	  int a;
	  double *weight = new double [2];
	  irtkMatrix *M = new irtkMatrix[2];

	  // check if time t fits a matrix
	  if ((int)round(t) == t) {
//		cout<<"using matrix "<<(int)round(t)<<endl;
//		for (int i = 0; i < 4; i++) {
//			for (int j = 0; j < 4; j++) {
//				cout<<_matrix[(int)round(t)](i, j)<<" , ";
//			}
//			cout<<endl;
//		}
		_tCurrent = time;
		_transformMatrix = _matrix[(int)round(t)];

		delete[] M;
		M = NULL;
		delete weight;
		weight = NULL;

		return _transformMatrix;
	  }

	  // compute weights for interpolation
	  a = (int)floor(t);
	  weight[1] = t - a;
	  weight[0] = 1 - weight[1];

	  // get matrizes for interpolation
	  M[0] = _matrix[a];
	  M[1] = _matrix[a+1];

	  ///////////////////////////////
//	  for (int x = 0; x < 4; x++) {
//	  cout<<_matrix[a](x, 0);
//	  for (int y = 1; y < 4; y++) {
//	  cout<<" , "<<_matrix[a](x, y);
//	  }
//	  cout<<endl;
//	  }
//	  for (int x = 0; x < 4; x++) {
//	  cout<<_matrix[a+1](x, 0);
//	  for (int y = 1; y < 4; y++) {
//	  cout<<" , "<<_matrix[a+1](x, y);
//	  }
//	  cout<<endl;
//	  }
	  ///////////////////////////////

	  // compute new matrix
	  _tCurrent = time;
	  _transformMatrix = FrechetMean(M, weight, 2);

	  delete[] M;
	  M = NULL;
	  delete weight;
	  weight = NULL;
	} else {
	  cerr << "wrong mode " << _mode << " -> interpolation not possible" << endl;
	  exit(1);
	}
  }

  // Set current time (and thus current matrix)
  _tCurrent = time;

  return _transformMatrix;
}

double irtkTemporalHomogeneousTransformation::Get(int index) const
{
  int i, j, t, indexT;

  if (index < this->NumberOfDOFs()) {
	t = index/12;
	indexT = index%12;
    i = indexT/4;
    j = indexT%4;
    return _matrix[t](i, j);
  } else {
    cerr << "irtkTemporalHomogeneousTransformation::Get: No such dof" << endl;
    return 0;
  }
}

double irtkTemporalHomogeneousTransformation::Get(int index, int t) const
{
  int i, j;

  if (index < this->NumberOfParams() && t < this->NumberOfTimePoints()) {
    i = index/4;
    j = index%4;
    return _matrix[t](i, j);
  } else {
    cerr << "irtkTemporalHomogeneousTransformation::Get: No such dof" << endl;
    return 0;
  }
}

double irtkTemporalHomogeneousTransformation::Get(int index, double time)
{
  int i, j;
  double t;

  t = this->TimeToLattice(time);
  this->GetMatrix(t);

  if (index < this->NumberOfParams() && t < this->NumberOfTimePoints()) {
    i = index/4;
    j = index%4;
    return _transformMatrix(i, j);
  } else {
    cerr << "irtkTemporalHomogeneousTransformation::Get: No such dof" << endl;
    return 0;
  }
}

void irtkTemporalHomogeneousTransformation::Put(int index, double x)
{
  int i, j, t, indexT;

  if (index < this->NumberOfDOFs()) {
	t = index/12;
	indexT = index%12;
	i = indexT/4;
	j = indexT%4;
    _matrix[t](i, j) = x;
  } else {
    cerr << "irtkTemporalHomogeneousTransformation::Put: No such dof" << endl;
    exit(1);
  }
}

void irtkTemporalHomogeneousTransformation::Put(int index, int t, double x)
{
  int i, j;

  if (index < this->NumberOfParams() && t < this->NumberOfTimePoints()) {
	i = index/4;
	j = index%4;
    _matrix[t](i, j) = x;
  } else {
    cerr << "irtkTemporalHomogeneousTransformation::Put: No such dof" << endl;
    exit(1);
  }
}

void irtkTemporalHomogeneousTransformation::Import(char *name)
{
  // Open file
  ifstream from(name, ios::in | ios::binary);

  // Check whether file opened ok
  if (!from) {
    cerr << "irtkTemporalHomogeneousTransformation::Import: Can't open file "
         << name << "\n";
    exit(1);
  }

  this->Import(from);
}

void irtkTemporalHomogeneousTransformation::Export(char *name)
{
  // Open file
  ofstream to(name, ios::out | ios::binary);

  // Check whether file opened ok
  if (!to) {
    cerr << "irtkTemporalHomogeneousTransformation::Export: Can't open file "
         << name << "\n";
    exit(1);
  }

  this->Export(to);
}

void irtkTemporalHomogeneousTransformation::Print()
{
  for (int i = 0; i < _Nt; i++) {
	cout<<"t: "<<i<<endl;
    _matrix[i].Print();
  }
}

void irtkTemporalHomogeneousTransformation::Transform(double &x, double &y, double &z, double t)
{
  double a, b, c;

  this->GetMatrix(t);

  // Pre-multiply point with transformation matrix
  a = (_transformMatrix(0, 0)*x+_transformMatrix(0, 1)*y+_transformMatrix(0, 2)*z+_transformMatrix(0, 3));
  b = (_transformMatrix(1, 0)*x+_transformMatrix(1, 1)*y+_transformMatrix(1, 2)*z+_transformMatrix(1, 3));
  c = (_transformMatrix(2, 0)*x+_transformMatrix(2, 1)*y+_transformMatrix(2, 2)*z+_transformMatrix(2, 3));

  // Copy result back
  x = a;
  y = b;
  z = c;
}

void irtkTemporalHomogeneousTransformation::Displacement(double &x, double &y, double &z, double t)
{
  double a, b, c;

  this->GetMatrix(t);
//  cout<<"Displacement matrix"<<endl;
//  cout<<M(0, 0)<<" , "<<M(0, 1)<<" , "<<M(0, 2)<<" , "<<M(0, 3)<<endl;
//  cout<<M(1, 0)<<" , "<<M(1, 1)<<" , "<<M(1, 2)<<" , "<<M(1, 3)<<endl;
//  cout<<M(2, 0)<<" , "<<M(2, 1)<<" , "<<M(2, 2)<<" , "<<M(2, 3)<<endl;

  // Pre-multiply point with transformation matrix
  a = (_transformMatrix(0, 0)*x+_transformMatrix(0, 1)*y+_transformMatrix(0, 2)*z+_transformMatrix(0, 3));
  b = (_transformMatrix(1, 0)*x+_transformMatrix(1, 1)*y+_transformMatrix(1, 2)*z+_transformMatrix(1, 3));
  c = (_transformMatrix(2, 0)*x+_transformMatrix(2, 1)*y+_transformMatrix(2, 2)*z+_transformMatrix(2, 3));

  // Copy result back
  x = a - x;
  y = b - y;
  z = c - z;
}

double irtkTemporalHomogeneousTransformation::Inverse(double &x, double &y, double &z, double t, double)
{
  double a, b, c;

  // Invert matrix
  irtkMatrix M = this->GetMatrix(t);
  M.Invert();

  // Pre-multiply point with transformation matrix
  a = (M(0, 0)*x+M(0, 1)*y+M(0, 2)*z+M(0, 3));
  b = (M(1, 0)*x+M(1, 1)*y+M(1, 2)*z+M(1, 3));
  c = (M(2, 0)*x+M(2, 1)*y+M(2, 2)*z+M(2, 3));

  // Copy result back
  x = a;
  y = b;
  z = c;

  // Error is zero by definition
  return 0;
}

bool irtkTemporalHomogeneousTransformation::IsIdentity()
{
  int i, j, t;

  for (t = 0; t < _Nt; t++) {
    for (i = 0; i < _matrix[t].Rows(); i++) {
      for (j = 0; j < _matrix[t].Cols(); j++) {
        if (i == j) {
          if (_matrix[t](i, j) != 1) return false;
        } else {
          if (_matrix[t](i, j) != 0) return false;
        }
      }
    }
  }
  return true;
}

bool irtkTemporalHomogeneousTransformation::IsIdentity(double t)
{
  int i, j;

  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  this->GetMatrix(t);

  for (i = 0; i < _transformMatrix.Rows(); i++) {
    for (j = 0; j < _transformMatrix.Cols(); j++) {
      if (i == j) {
        if (_transformMatrix(i, j) != 1) return false;
      } else {
        if (_transformMatrix(i, j) != 0) return false;
      }
    }
  }
  return true;
}

void irtkTemporalHomogeneousTransformation::Invert()
{
  // Invert transformations
  for (int t = 0; t < _Nt; t++)
    _matrix[t].Invert();
}

void irtkTemporalHomogeneousTransformation::Invert(int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  // Invert transformation
  _matrix[t].Invert();
}

void irtkTemporalHomogeneousTransformation::Jacobian(irtkMatrix &jac, double, double, double, double t)
{
  int i, j;

  this->GetMatrix(t);

  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Calculate matrix
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      jac(i, j) = _transformMatrix(i, j);
    }
  }
}

void irtkTemporalHomogeneousTransformation::GlobalJacobian(irtkMatrix &jac, double, double, double, double t)
{
  int i, j;

  this->GetMatrix(t);

  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Calculate matrix
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      jac(i, j) = _transformMatrix(i, j);
    }
  }
}

void irtkTemporalHomogeneousTransformation::LocalJacobian(irtkMatrix &jac, double, double, double, double t)
{
  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Set matrix to identity
  jac(0, 0) = 1;
  jac(1, 1) = 1;
  jac(2, 2) = 1;
}

void irtkTemporalHomogeneousTransformation::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double t)
{
  int i, j;

  // calculate matrix indizes
  i = dof/4;
  j = dof%4;

  jac[0] = jac[1] = jac[2] = 0.;

  switch (i) {
  case 0:
	switch (j) {
	case 0:
	  jac[0] = x;
	  break;
	case 1:
	  jac[0] = y;
	  break;
	case 2:
	  jac[0] = z;
	  break;
	case 3:
	  jac[0] = 1.;
	  break;
	default:
	  cerr << this->NameOfClass() << "::JacobianDOFs(): No such dof = " << dof << endl;
	  exit(1);
	  break;
	}
	break;
  case 1:
	switch (j) {
	case 0:
	  jac[1] = x;
	  break;
	case 1:
	  jac[1] = y;
	  break;
	case 2:
	  jac[1] = z;
	  break;
	case 3:
	  jac[1] = 1.;
	  break;
	default:
	  cerr << this->NameOfClass() << "::JacobianDOFs(): No such dof = " << dof << endl;
	  exit(1);
	  break;
	}
	break;
  case 2:
	switch (j) {
	case 0:
	  jac[2] = x;
	  break;
	case 1:
	  jac[2] = y;
	  break;
	case 2:
	  jac[2] = z;
	  break;
	case 3:
	  jac[2] = 1.;
	  break;
	default:
	  cerr << this->NameOfClass() << "::JacobianDOFs(): No such dof = " << dof << endl;
	  exit(1);
	  break;
	}
	break;
  default:
	cerr << this->NameOfClass() << "::JacobianDOFs(): No such dof = " << dof << endl;
	exit(1);
	break;
  }
}

void irtkTemporalHomogeneousTransformation::SubdivideOdd()
{
  int k;
  double i, t;

  cerr << "matrix spacing might not be frame spacing after this -> optimisation impossible" << endl;
  cerr << "check before and use SubdivideEven instead" << endl;

  // Allocate memory for new control points
  irtkMatrix * matrix = NULL;
  matrix = this->Allocate(matrix, 2*_Nt-1);

  k = 0;
  for (i = 0; i < _Nt-0.5; i+=0.5) {
	t = LatticeToTime(i);
    matrix[k] = GetMatrix(t);
    k++;
  }

  // Deallocate points
  _matrix = this->Deallocate(_matrix, _Nt);
  delete []_status;

  // Update pointers to control points
  _matrix = matrix;

  // Increase number of control points
  _Nt = 2*_Nt-1;

  // Recalculate control point spacing
  _dt = _dt/2;

  // Initialize memory for control point status
  _status = new _Status[this->NumberOfDOFs()];
  for (k = 0; k < this->NumberOfDOFs(); k++) {
    _status[k] = _Active;
  }
}

void irtkTemporalHomogeneousTransformation::SubdivideEven()
{
  int k;
  double i, t;

  cerr << "matrix spacing might not be frame spacing after this -> optimisation impossible" << endl;
  cerr << "check before and use SubdivideOdd instead" << endl;

  // Allocate memory for new control points
  irtkMatrix * matrix = NULL;
  matrix = this->Allocate(matrix, 2*_Nt);

  k = 0;
  for (i = 0; i < _Nt; i+=0.5) {
	t = LatticeToTime(i);
    matrix[k] = GetMatrix(t);
    k++;
  }

  // Deallocate points
  _matrix = this->Deallocate(_matrix, _Nt);
  delete []_status;

  // Update pointers to control points
  _matrix = matrix;

  // Increase number of control points
  _Nt = 2*_Nt;

  // Recalculate control point spacing
  _dt = _dt/3;

  // Initialize memory for control point status
  _status = new _Status[this->NumberOfDOFs()];
  for (k = 0; k < this->NumberOfDOFs(); k++) {
    _status[k] = _Active;
  }
}

istream& irtkTemporalHomogeneousTransformation::Import(istream& is)
{
  for (int t = 0; t < _Nt; t++) {
    is >> this->_matrix[t];

    // Invert transformation matrix
    this->Invert(t);
  }

  return is;
}

ostream& irtkTemporalHomogeneousTransformation::Export(ostream& os)
{
  for (int t = 0; t < _Nt; t++) {
    // Invert transformation matrix
    this->Invert(t);

    os << this->_matrix[t];

    // Invert transformation matrix
    this->Invert(t);
  }

  return os;
}

irtkCifstream& irtkTemporalHomogeneousTransformation::Read(irtkCifstream& from)
{
  int i;
  double data, dt, tMin, tMax;
  unsigned int magic_no, trans_type, dofs, Nt;

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkTemporalHomogeneousTransformation::Read: Not a vaild transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if (trans_type != IRTKTRANSFORMATION_HOMO_TEMPORAL) {
    cerr << "irtkTemporalHomogeneousTransformation::Read: Not a vaild affine transformation" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&dofs, 1);
  if (dofs != 12) {
    cerr << "irtkTemporalHomogeneousTransformation::Read: Invalid no. of dofs = " << dofs << endl;
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

irtkCofstream& irtkTemporalHomogeneousTransformation::Write(irtkCofstream& to)
{
  int i;

  // Write magic no. for transformations
  unsigned int magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  unsigned int trans_type = IRTKTRANSFORMATION_HOMO_TEMPORAL;
  to.WriteAsUInt(&trans_type, 1);

  // Write transformation type
  unsigned int dofs = this->NumberOfDOFs()/this->_Nt;
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

double *irtkTemporalHomogeneousTransformation::Allocate(double *data, int t)
{
  int i;

  if (t == 0) {
    return NULL;
  }

  if ((data = new double [t+8]) == NULL) {
    cerr << "Allocate: malloc failed for " << t << "\n";
    exit(1);
  }
  data += 4;

  for (i = -4; i < t+4; i++) {
	data[i] = 0;
  }

  return data;
}

double *irtkTemporalHomogeneousTransformation::Deallocate(double *data, int)
{
  if (data != NULL) {
    delete [](data-4);
  }
  return NULL;
}

irtkMatrix *irtkTemporalHomogeneousTransformation::Allocate(irtkMatrix *data, int t)
{
  int i;

  if (t == 0) {
    return NULL;
  }

  if ((data = new irtkMatrix [t+8]) == NULL) {
    cerr << "Allocate: malloc failed for " << t << "\n";
    exit(1);
  }
  data += 4;

  for (i = -4; i < t+4; i++) {
	data[i] = irtkMatrix(4, 4);
	data[i].Ident();
  }

  return data;
}

irtkMatrix *irtkTemporalHomogeneousTransformation::Deallocate(irtkMatrix *data, int)
{
  if (data != NULL) {
    delete [](data-4);
  }
  return NULL;
}




