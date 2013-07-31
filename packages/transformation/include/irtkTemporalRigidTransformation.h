
#ifndef _IRTKTEMPORALRIGIDTRANSFORMATION_H

#define _IRTKTEMPORALRIGIDTRANSFORMATION_H

/**
 * Class for rigid transformations.
 *
 * This class defines and implements rigid body transformations. The rigid
 * body transformations are parameterized by three rotations around the axes
 * of the coordinate system followed by three translations along the axes of
 * the coordinate system. Note that the order of rotations is defined as a
 * rotation around the z-axis, the y-axis and finally around the x-axis. In
 * total, the transformation is parameterized by six degrees of freedom.
 *
 */

class irtkTemporalRigidTransformation : public irtkTemporalHomogeneousTransformation
{

protected:

  /// Translation along the x-axis (in mm)
  double *_tx;

  /// Translation along the y-axis (in mm)
  double *_ty;

  /// Translation along the z-axis (in mm)
  double *_tz;

  /// Rotation around the x-axis (in degrees)
  double *_rx;

  /// Rotation around the y-axis (in degrees)
  double *_ry;

  /// Rotation around the z-axis (in degrees)
  double *_rz;

  /// Construct a matrix based on parameters passed in the array.
  virtual irtkMatrix Parameters2Matrix(double *) const;

  /// Return an array with parameters corresponding to a given matrix.
  virtual void Matrix2Parameters(irtkMatrix, double*) const;

  /// Updates transformation matrizes (at control points)
  virtual void UpdateMatrix();
  virtual void UpdateMatrix(int);

public:

  /// Constructor (default)
  irtkTemporalRigidTransformation();
  irtkTemporalRigidTransformation(double, double, double);

  /// Constructor (copy)
  irtkTemporalRigidTransformation(const irtkTemporalRigidTransformation &);

  /// Destructor
  virtual ~irtkTemporalRigidTransformation();

  virtual void Initialize(double, double, double);

  /// Reset transformation
  virtual void Reset();
  virtual void Reset(int);

  /// Puts translation along the x-axis (control point)
  void   PutTranslationX(double, int);

  /// Gets translation along the x-axis
  virtual double GetTranslationX(int) const;

  /// Puts translation along the y-axis (control point)
  virtual void   PutTranslationY(double, int);

  /// Gets translation along the y-axis
  virtual double GetTranslationY(int) const;

  /// Puts translation along the z-axis (control point)
  virtual void   PutTranslationZ(double, int);

  /// Gets translation along the z-axis
  virtual double GetTranslationZ(int) const;

  /// Puts rotation angle around the x-axis (control point)
  virtual void   PutRotationX(double, int);

  /// Gets rotation angle around the x-axis
  virtual double GetRotationX(int) const;

  /// Puts rotation angle around the y-axis (control point)
  virtual void   PutRotationY(double, int);

  /// Gets rotation angle around the y-axis
  virtual double GetRotationY(int) const;

  /// Puts rotation angle around the z-axis (control point)
  virtual void   PutRotationZ(double, int);

  /// Gets rotation angle around the z-axis
  virtual double GetRotationZ(int) const;

  /// Returns the number of degrees of freedom of the transformation
  virtual int NumberOfDOFs() const;

  /// Returns the number of parameters of the transformation
  virtual int NumberOfParams() const;

  /// Puts the transformation matrix (control points) (Argument is not checked)
  virtual void PutMatrix(const irtkMatrix &, int);

  /// Puts a transformation parameter (control point)
  virtual void   Put(int, double);
  virtual void   Put(int, int, double);

  /// Gets a transformation parameter (control point)
  virtual double Get(int) const;
  virtual double Get(int, int) const;
  virtual double Get(int, double);

  /// Assign the parameters passed to the current object (control point values)
  virtual void SetParameters(double **params);
  virtual void SetParameters(double *params, int);

  /// Transforms a point by the rotation part of the rigid transformation.
  virtual void Rotate(double& x, double& y, double& z, double = 0);

  /// Calculate the Jacobian of the transformation with respect to the transformation parameters
  /// B-splines only in the temporal direction, so no after-differentiation (right???)
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0);

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity();
  virtual bool IsIdentity(int);

  /// Subdivide temporal control point spacing (for odd numbers of frames)
  virtual void SubdivideOdd();

  /// Subdivide temporal control point spacing (for even numbers of frames)
  virtual void SubdivideEven();

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Check file header
  static int CheckHeader(char *);

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

  /// Reads a transformation from a file
  virtual irtkCifstream& Read(irtkCifstream&);

  /// Writes a transformation to a file
  virtual irtkCofstream& Write(irtkCofstream&);

  /// Imports a transformation from a file
  virtual istream& Import(istream&);

  /// Exports a transformation to a file
  virtual ostream& Export(ostream&);
};

inline int irtkTemporalRigidTransformation::NumberOfDOFs() const
{
  return _Nt*6;
}

inline int irtkTemporalRigidTransformation::NumberOfParams() const
{
  return 6;
}

inline irtkTemporalRigidTransformation::irtkTemporalRigidTransformation()
{
  int i;

  // Initialize control point parameters (time)
  _Nt = 2;
  _tMin = 0;
  _tMax = 1;
  _dt = 1;

  _mode = Linear;

  _tx = this->Allocate(_tx, _Nt);
  _ty = this->Allocate(_ty, _Nt);
  _tz = this->Allocate(_tz, _Nt);
  _rx = this->Allocate(_rx, _Nt);
  _ry = this->Allocate(_ry, _Nt);
  _rz = this->Allocate(_rz, _Nt);

  _matrix = this->Allocate(_matrix, _Nt);

  // Free memory allocated for DOF status by base class
  delete []_status;

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = _Active;
  }

  _transformMatrix = _matrix[0];
  _tCurrent = _tMin;
}

inline irtkTemporalRigidTransformation::irtkTemporalRigidTransformation(double dt, double t1, double t2)
{
  int i;

  // Initialize control point parameters (time)
  _tMin = t1;
  _tMax = t2;
  _Nt = round((_tMax - _tMin) / dt) + 1;
  _dt = (_tMax - _tMin) / (_Nt - 1);

  _mode = Linear;

  _tx = this->Allocate(_tx, _Nt);
  _ty = this->Allocate(_ty, _Nt);
  _tz = this->Allocate(_tz, _Nt);
  _rx = this->Allocate(_rx, _Nt);
  _ry = this->Allocate(_ry, _Nt);
  _rz = this->Allocate(_rz, _Nt);

  _matrix = this->Allocate(_matrix, _Nt);

  // Free memory allocated for DOF status by base class
  delete []_status;

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = _Active;
  }

  _transformMatrix = _matrix[0];
  _tCurrent = _tMin;
}

inline irtkTemporalRigidTransformation::irtkTemporalRigidTransformation(const irtkTemporalRigidTransformation &t) : irtkTemporalHomogeneousTransformation(t)
{
  int i;

  _tMin = t.GetTMin();
  _tMax = t.GetTMax();
  _Nt = t.NumberOfTimePoints();
  _dt = t.GetDT();

  _mode = Linear;

  _tx = this->Allocate(_tx, _Nt);
  _ty = this->Allocate(_ty, _Nt);
  _tz = this->Allocate(_tz, _Nt);
  _rx = this->Allocate(_rx, _Nt);
  _ry = this->Allocate(_ry, _Nt);
  _rz = this->Allocate(_rz, _Nt);

  for (i = -4; i < _Nt+4; i++) {
    _tx[i] = t._tx[i];
    _ty[i] = t._ty[i];
    _tz[i] = t._tz[i];
    _rx[i] = t._rx[i];
    _ry[i] = t._ry[i];
    _rz[i] = t._rz[i];
  }

  _matrix = this->Allocate(_matrix, _Nt);

  // Free memory allocated for DOF status by base class
  delete []_status;

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = t._status[i];
  }

  this->UpdateMatrix();
  _transformMatrix = _matrix[0];
  _tCurrent = _tMin;
}

inline irtkTemporalRigidTransformation::~irtkTemporalRigidTransformation()
{
  // Free memory for control points if necessary
  if (_tx != NULL) _tx = this->Deallocate(_tx, _Nt);
  if (_ty != NULL) _ty = this->Deallocate(_ty, _Nt);
  if (_tz != NULL) _tz = this->Deallocate(_tz, _Nt);
  if (_rx != NULL) _rx = this->Deallocate(_rx, _Nt);
  if (_ry != NULL) _ry = this->Deallocate(_ry, _Nt);
  if (_rz != NULL) _rz = this->Deallocate(_rz, _Nt);

  if (_matrix != NULL) _matrix = this->Deallocate(_matrix, _Nt);

  _Nt = 0;
}

inline void irtkTemporalRigidTransformation::Initialize(double dt, double t1, double t2)
{
  int i;

  // Free memory for control points if necessary
  if (_tx != NULL) _tx = this->Deallocate(_tx, _Nt);
  if (_ty != NULL) _ty = this->Deallocate(_ty, _Nt);
  if (_tz != NULL) _tz = this->Deallocate(_tz, _Nt);
  if (_rx != NULL) _rx = this->Deallocate(_rx, _Nt);
  if (_ry != NULL) _ry = this->Deallocate(_ry, _Nt);
  if (_rz != NULL) _rz = this->Deallocate(_rz, _Nt);

  if (_matrix != NULL) _matrix = this->Deallocate(_matrix, _Nt);

  // Initialize control point parameters (time)
  _tMin = t1;
  _tMax = t2;
  _Nt = round((_tMax - _tMin) / dt) + 1;
  _dt = (_tMax - _tMin) / (_Nt - 1);

  _mode = Linear;

  _tx = this->Allocate(_tx, _Nt);
  _ty = this->Allocate(_ty, _Nt);
  _tz = this->Allocate(_tz, _Nt);
  _rx = this->Allocate(_rx, _Nt);
  _ry = this->Allocate(_ry, _Nt);
  _rz = this->Allocate(_rz, _Nt);

  _matrix = this->Allocate(_matrix, _Nt);

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = _Active;
  }

  _transformMatrix = _matrix[0];
  _tCurrent = _tMin;
}

inline void irtkTemporalRigidTransformation::Reset()
{
  for (int i = 0; i < _Nt; i++) {
    this->Reset(i);
  }
}

inline void irtkTemporalRigidTransformation::Reset(int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  // Initialize rotations and translations
  _tx[t] = _ty[t] = _tz[t] = 0;
  _rx[t] = _ry[t] = _rz[t] = 0;
  this->UpdateMatrix(t);
}

inline void irtkTemporalRigidTransformation::PutMatrix(const irtkMatrix &matrix, int i)
{
  double *params = new double [6];

  if (i < 0 || i >= _Nt)
	cerr << "wrong index "<<i<<endl;

  _matrix[i] = matrix;

  this->Matrix2Parameters(matrix, params);
  _tx[i] = params[TX];
  _ty[i] = params[TY];
  _tz[i] = params[TZ];
  _rx[i] = params[RX];
  _ry[i] = params[RY];
  _rz[i] = params[RZ];
}

inline void irtkTemporalRigidTransformation::PutRotationX(double rx, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _rx[t] = rx;
  this->UpdateMatrix(t);
}

inline void irtkTemporalRigidTransformation::PutRotationY(double ry, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _ry[t] = ry;
  this->UpdateMatrix(t);
}

inline void irtkTemporalRigidTransformation::PutRotationZ(double rz, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _rz[t] = rz;
  this->UpdateMatrix(t);
}

inline void irtkTemporalRigidTransformation::PutTranslationX(double tx, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _tx[t] = tx;
  this->UpdateMatrix(t);
}

inline void irtkTemporalRigidTransformation::PutTranslationY(double ty, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _ty[t] = ty;
  this->UpdateMatrix(t);
}

inline void irtkTemporalRigidTransformation::PutTranslationZ(double tz, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _tz[t] = tz;
  this->UpdateMatrix(t);
}

inline const char *irtkTemporalRigidTransformation::NameOfClass()
{
  return "irtkTemporalRigidTransformation";
}

inline void irtkTemporalRigidTransformation::SetParameters(double **params)
{
  for (int i = 0; i < _Nt; i++) {
    _tx[i] = params[i][TX];
    _ty[i] = params[i][TY];
    _tz[i] = params[i][TZ];

    _rx[i] = params[i][RX];
    _ry[i] = params[i][RY];
    _rz[i] = params[i][RZ];
  }
  this->UpdateMatrix();
}

inline void irtkTemporalRigidTransformation::SetParameters(double *params, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _tx[t] = params[TX];
  _ty[t] = params[TY];
  _tz[t] = params[TZ];

  _rx[t] = params[RX];
  _ry[t] = params[RY];
  _rz[t] = params[RZ];
  this->UpdateMatrix(t);
}

#endif
