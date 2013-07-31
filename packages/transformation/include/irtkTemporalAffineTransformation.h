
#ifndef _IRTKTEMPORALAFFINETRANSFORMATION_H

#define _IRTKTEMPORALAFFINETRANSFORMATION_H

/**
 * Class for affine transformations.
 *
 * This class defines and implements affine transformations. In addition to
 * the rigid body transformation parameters, affine transformations are
 * parameterized by three scaling and and six skewing parameters. The three
 * scaling parameters define the scaling along the axis of the coordinate
 * transformations. The six skewing parameters define the skewing angles in
 * different planes. Note that the six skewing parameters are not independent.
 * In total, the transformation can be parameterized by up to 15 degrees of
 * freedom of which only 12 are linearly independent. The class in its current
 * implementation supports 9, 12 or 15 degrees of freedom. The number of
 * degrees of freedom can be specified when constructing the transformation.
 * The default is 15 degrees of freedom.
 *
 */

class irtkTemporalAffineTransformation : public irtkTemporalRigidTransformation
{

protected:

  /// Scaling along the x-axis
  double *_sx;

  /// Scaling along the y-axis
  double *_sy;

  /// Scaling along the z-axis
  double *_sz;

  /// Skew angle in the x direction based on y component (in degrees)
  double *_sxy;

  /// Skew angle in the y direction based on z component (in degrees)
  double *_syz;

  /// Skew angle in the x direction based on z component (in degrees)
  double *_sxz;

  /// Construct a matrix based on parameters passed in the array.
  virtual irtkMatrix Parameters2Matrix(double *) const;

  /// Return an array with parameters corresponding to a given matrix.
  virtual void Matrix2Parameters(irtkMatrix, double*) const;

  /// Updates transformation matrizes (at control points)
  virtual void UpdateMatrix();
  virtual void UpdateMatrix(int);

public:

  /// Constructor (default)
  irtkTemporalAffineTransformation();
  irtkTemporalAffineTransformation(double, double, double);

  /// Constructor (copy)
  irtkTemporalAffineTransformation(const irtkTemporalRigidTransformation &);

  /// Constructor (copy)
  irtkTemporalAffineTransformation(const irtkTemporalAffineTransformation &);

  /// Destructor
  virtual ~irtkTemporalAffineTransformation();

  virtual void Initialize(double, double, double);

  /// Reset transformation
  virtual void Reset();
  virtual void Reset(int);

  /// Puts scaling factor along the x-axis (control point)
  virtual void   PutScaleX(double, int);

  /// Gets scaling factor along the x-axis
  virtual double GetScaleX(int) const;

  /// Puts scaling factor along the y-axis (control point)
  virtual void   PutScaleY(double, int);

  /// Gets scaling factor along the y-axis
  virtual double GetScaleY(int) const;

  /// Puts scaling factor along the z-axis (control point)
  virtual void   PutScaleZ(double, int);

  /// Gets scaling factor along the z-axis
  virtual double GetScaleZ(int) const;

  /// Puts y-dependent skewing angle in the x direction (in degrees) (control point)
  virtual void   PutShearXY(double, int);

  /// Gets y-dependent skewing angle in the x direction (in degrees)
  virtual double GetShearXY(int) const;

  /// Puts z-dependent skewing angle in the y direction (in degrees) (control point)
  virtual void   PutShearYZ(double, int);

  /// Gets z-dependent skewing angle in the y direction (in degrees)
  virtual double GetShearYZ(int) const;

  /// Puts z-dependent skewing angle in the x direction (in degrees) (control point)
  virtual void   PutShearXZ(double, int);

  /// Gets z-dependent skewing angle in the x direction (in degrees)
  virtual double GetShearXZ(int) const;

  /// Returns the number of degrees of freedom of the transformation
  virtual int NumberOfDOFs() const;

  /// Returns the number of parameters of the transformation
  virtual int NumberOfParams() const;

  /// Puts the transformation matrix (control points) (Argument is not checked)
  virtual void PutMatrix(const irtkMatrix &, int);

  /// Puts a transformation parameter (transformation matrix is updated)
  virtual void   Put(int, double);
  virtual void   Put(int, int, double);

  /// Gets a transformation parameter
  virtual double Get(int) const;
  virtual double Get(int, int) const;
  virtual double Get(int, double);

  /// Assign the parameters passed in the array to the current object (control point values)
  virtual void SetParameters(double **);
  virtual void SetParameters(double *, int);

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

inline int irtkTemporalAffineTransformation::NumberOfDOFs() const
{
  return _Nt*12;
}

inline int irtkTemporalAffineTransformation::NumberOfParams() const
{
  return 12;
}

inline irtkTemporalAffineTransformation::irtkTemporalAffineTransformation()
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
  _sx = this->Allocate(_sx, _Nt);
  _sy = this->Allocate(_sy, _Nt);
  _sz = this->Allocate(_sz, _Nt);
  _sxy = this->Allocate(_sxy, _Nt);
  _sxz = this->Allocate(_sxz, _Nt);
  _syz = this->Allocate(_syz, _Nt);

  for (i = -4; i < _Nt+4; i++) {
	_sx[i] = 100;
	_sy[i] = 100;
	_sz[i] = 100;
  }

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

inline irtkTemporalAffineTransformation::irtkTemporalAffineTransformation(double dt, double t1, double t2)
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
  _sx = this->Allocate(_sx, _Nt);
  _sy = this->Allocate(_sy, _Nt);
  _sz = this->Allocate(_sz, _Nt);
  _sxy = this->Allocate(_sxy, _Nt);
  _sxz = this->Allocate(_sxz, _Nt);
  _syz = this->Allocate(_syz, _Nt);

  for (i = -4; i < _Nt+4; i++) {
	_sx[i] = 100;
	_sy[i] = 100;
	_sz[i] = 100;
  }

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

inline irtkTemporalAffineTransformation::irtkTemporalAffineTransformation(const irtkTemporalRigidTransformation &t)
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
  _sx = this->Allocate(_sx, _Nt);
  _sy = this->Allocate(_sy, _Nt);
  _sz = this->Allocate(_sz, _Nt);
  _sxy = this->Allocate(_sxy, _Nt);
  _sxz = this->Allocate(_sxz, _Nt);
  _syz = this->Allocate(_syz, _Nt);

  for (i = -4; i < _Nt+4; i++) {
    // Copy rotations and translations
    _tx[i] = t.Get(0, i);
    _ty[i] = t.Get(1, i);
    _tz[i] = t.Get(2, i);
    _rx[i] = t.Get(3, i);
    _ry[i] = t.Get(4, i);
    _rz[i] = t.Get(5, i);

	// Initialize scale and shears
    _sx[i]  = 100;
    _sy[i]  = 100;
    _sz[i]  = 100;
  }

  _matrix = this->Allocate(_matrix, _Nt);

  // Free memory allocated for DOF status by base class
  delete []_status;

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    if (i < 6) {
      _status[i] = t.GetStatus(i);
    } else {
      _status[i] = _Active;
    }
  }

  this->UpdateMatrix();
  _transformMatrix = _matrix[0];
  _tCurrent = _tMin;
}

inline irtkTemporalAffineTransformation::irtkTemporalAffineTransformation(const irtkTemporalAffineTransformation &t) : irtkTemporalRigidTransformation(t)
{
  int i;

  _tMin = t._tMin;
  _tMax = t._tMax;
  _Nt = t._Nt;
  _dt = t._dt;

  _mode = Linear;

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

  for (i = -4; i < _Nt+4; i++) {
    // Copy rotations and translations
	_tx[i] = t._tx[i];
	_ty[i] = t._ty[i];
	_tz[i] = t._tz[i];
	_rx[i] = t._rx[i];
	_ry[i] = t._ry[i];
	_rz[i] = t._rz[i];

	// Copy scale and shears
	_sx[i]  = t._sx[i];
	_sy[i]  = t._sy[i];
	_sz[i]  = t._sz[i];
	_sxy[i] = t._sxy[i];
	_syz[i] = t._syz[i];
	_sxz[i] = t._sxz[i];
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

inline irtkTemporalAffineTransformation::~irtkTemporalAffineTransformation()
{
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

  if (_matrix != NULL) _matrix = this->Deallocate(_matrix, _Nt);

  _Nt = 0;
}

inline void irtkTemporalAffineTransformation::Initialize(double dt, double t1, double t2)
{
  int i;

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
  _sx = this->Allocate(_sx, _Nt);
  _sy = this->Allocate(_sy, _Nt);
  _sz = this->Allocate(_sz, _Nt);
  _sxy = this->Allocate(_sxy, _Nt);
  _sxz = this->Allocate(_sxz, _Nt);
  _syz = this->Allocate(_syz, _Nt);

  for (i = -4; i < _Nt+4; i++) {
	_sx[i] = 100;
	_sy[i] = 100;
	_sz[i] = 100;
  }

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

inline void irtkTemporalAffineTransformation::Reset()
{
  for (int i = 0; i < _Nt; i++) {
	this->Reset(i);
  }
}

inline void irtkTemporalAffineTransformation::Reset(int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  // Reset rigid part
  this->irtkTemporalRigidTransformation::Reset(t);

  // Initialize scale and shears
  _sx[t]  = _sy[t]  = _sz[t]  = 100;
  _sxy[t] = _syz[t] = _sxz[t] = 0;
  this->UpdateMatrix(t);
}

inline void irtkTemporalAffineTransformation::PutMatrix(const irtkMatrix &matrix, int i)
{
  double *params = new double [12];

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
  _sx[i] = params[SX];
  _sy[i] = params[SY];
  _sz[i] = params[SZ];
  _sxy[i] = params[SXY];
  _sxz[i] = params[SXZ];
  _syz[i] = params[SYZ];
}

inline void irtkTemporalAffineTransformation::PutScaleX(double sx, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _sx[t] = sx;
  this->UpdateMatrix(t);
}

inline void irtkTemporalAffineTransformation::PutScaleY(double sy, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _sy[t] = sy;
  this->UpdateMatrix(t);
}

inline void irtkTemporalAffineTransformation::PutScaleZ(double sz, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _sz[t] = sz;
  this->UpdateMatrix(t);
}

inline void irtkTemporalAffineTransformation::PutShearXY(double sxy, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _sxy[t] = sxy;
  this->UpdateMatrix(t);
}

inline void irtkTemporalAffineTransformation::PutShearYZ(double syz, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _syz[t] = syz;
  this->UpdateMatrix(t);
}

inline void irtkTemporalAffineTransformation::PutShearXZ(double sxz, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  _sxz[t] = sxz;
  this->UpdateMatrix(t);
}

inline const char *irtkTemporalAffineTransformation::NameOfClass()
{
  return "irtkTemporalAffineTransformation";
}

inline void irtkTemporalAffineTransformation::SetParameters(double **params)
{
  this->irtkTemporalRigidTransformation::SetParameters(params);

  for (int i = 0; i < _Nt; i++) {
    _sx[i]  = params[i][SX];
    _sy[i]  = params[i][SY];
    _sz[i]  = params[i][SZ];

    _sxy[i] = params[i][SXY];
    _syz[i] = params[i][SYZ];
    _sxz[i] = params[i][SXZ];
  }
  this->UpdateMatrix();
}

inline void irtkTemporalAffineTransformation::SetParameters(double *params, int t)
{
  if (t < 0 || t >= _Nt)
	cerr << "wrong index "<<t<<endl;

  this->irtkTemporalRigidTransformation::SetParameters(params, t);

  _sx[t]  = params[SX];
  _sy[t]  = params[SY];
  _sz[t]  = params[SZ];

  _sxy[t] = params[SXY];
  _syz[t] = params[SYZ];
  _sxz[t] = params[SXZ];
  this->UpdateMatrix(t);
}

#endif
