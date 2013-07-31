

#ifndef _TEMPORALHOMOGENEOUSTRANSFORMATION_H

#define _TEMPORALHOMOGENEOUSTRANSFORMATION_H

/**
 * Base class for homogeneous transformations.
 *
 * This class defines and implements homogeneous transformations which can
 * be represented by a 4 x 4 transformation matrix. The transformation of
 * a point is implemented by post multiplying the transformation matrix with
 * the point in homogenous coordinates. The transformation is parameterized
 * by twelve degrees of freedom.
 *
 */

//typedef enum {TX, TY, TZ, RX, RY, RZ, SX, SY, SZ, SXY, SYZ, SXZ, SYX, SZY, SZX}
//irtkTemporalBSplineHomogeneousTransformationParameterIndex;

typedef enum {Linear, CSpline}
irtkTemporalHomogeneousTransformationInterpolationIndex;

class irtkTemporalHomogeneousTransformation : public irtkTransformation
{

protected:

  /// matrix interpolation mode
  int _mode;

  /// 4 x 4 transformation matrix for homogeneous coordinates
  irtkMatrix *_matrix;

  /// 4 x 4 transformation matrix at time t (used to safe interpolated matrix -> efficiency)
  irtkMatrix _transformMatrix;

  /// time point of current _transformation matrix
  double _tCurrent;

  /// Number of temporal control points
  int _Nt;

  /// Spacing of control points in t (in ms)
  double _dt;

  /// The minimum time value for which the transformation is defined.
  double _tMin;

  /// The maximum time value for which the transformation is defined.
  double _tMax;

  /// Allocate memory for control points
  static double *Allocate  (double *, int);

  /// Deallocate memory for control points
  static double *Deallocate(double *, int);

  /// Allocate memory for control points
  static irtkMatrix *Allocate  (irtkMatrix *, int);

  /// Deallocate memory for control points
  static irtkMatrix *Deallocate(irtkMatrix *, int);

public:

  /// Constructor (default)
  irtkTemporalHomogeneousTransformation();
  irtkTemporalHomogeneousTransformation(double, double, double);

  /// Constructor (from matrix)
  irtkTemporalHomogeneousTransformation(const irtkMatrix *, double, double, double);

  /// Constructor (copy)
  irtkTemporalHomogeneousTransformation(const irtkTemporalHomogeneousTransformation &);

  /// Destructor
  virtual ~irtkTemporalHomogeneousTransformation();

  virtual void Initialize(double, double, double);

  /// Returns the number of degrees of freedom of the transformation
  virtual int    NumberOfDOFs() const;

  /// Returns the number of time points of the transformation
  virtual int    NumberOfTimePoints() const;

  /// Returns the number of parameters of the transformation
  virtual int    NumberOfParams() const;

  /// Returns delta T
  virtual double GetDT() const;

  /// Returns minimum of time interval
  virtual double GetTMin() const;

  /// Returns maximum of time interval
  virtual double GetTMax() const;

  /// Set interpolation mode
  virtual void SetIntMode(int);

  /// Get interpolation mode
  virtual int GetInMode() const;

  /// Gets a transformation parameter (control point)
  virtual double Get(int) const;
  virtual double Get(int, int) const;
  virtual double Get(int, double);

  /// Puts a transformation paramater (control point)
  virtual void   Put(int, double);
  virtual void   Put(int, int, double);

  /// Gets transformation matrix at time t
  virtual irtkMatrix GetMatrix(double);

  /// Puts the transformation matrix (control points) (Argument is not checked)
  virtual void PutMatrix(const irtkMatrix &, int);
  virtual void PutMatrix(irtkMatrix *);

  /// Transforms a single point
  virtual void Transform(double &, double &, double &, double = 0);

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0);

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double &, double &, double &, double = 0);

  /// Calculates displacement
  virtual void Displacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &, double = 0);

  /// Inverts the transformation matrix at time t
  virtual void Invert();
  virtual void Invert(int);

  /// Inverse transformation
  virtual double Inverse(double &, double &, double &, double = 0, double = 0.01);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the global transformation
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0);

  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0);

  /// Checks whether transformation is an identity mapping
  virtual bool IsIdentity();
  virtual bool IsIdentity(double t);

  /// Subdivide temporal control point spacing (for odd numbers of frames)
  virtual void SubdivideOdd();

  /// Subdivide temporal control point spacing (for even numbers of frames)
  virtual void SubdivideEven();

  /// Lattice to time transformation
  virtual double LatticeToTime(double) const;

  /// Time to lattice transformation
  virtual double TimeToLattice(double) const;

  /// Prints the parameters of the transformation
  virtual void Print();

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

  /// Reads a transformation from a file
  virtual irtkCifstream& Read(irtkCifstream&);

  /// Writes a transformation to a file
  virtual irtkCofstream& Write(irtkCofstream&);

  /// Imports a transformation from a file
  virtual void Import(char *);

  /// Exports a transformation to a file
  virtual void Export(char *);

  /// Imports a transformation from a file
  virtual istream& Import(istream&);

  /// Exports a transformation to a file
  virtual ostream& Export(ostream&);

};

inline int irtkTemporalHomogeneousTransformation::NumberOfDOFs() const
{
  return 12*_Nt;
}

inline int irtkTemporalHomogeneousTransformation::NumberOfTimePoints() const
{
  return _Nt;
}

inline int irtkTemporalHomogeneousTransformation::NumberOfParams() const
{
  return 12;
}

inline double irtkTemporalHomogeneousTransformation::GetDT() const
{
  return _dt;
}

inline double irtkTemporalHomogeneousTransformation::GetTMin() const
{
  return _tMin;
}

inline double irtkTemporalHomogeneousTransformation::GetTMax() const
{
  return _tMax;
}

inline irtkTemporalHomogeneousTransformation::irtkTemporalHomogeneousTransformation()
{
  int i;

  // Initialize control point parameters (time)
  _Nt = 2;
  _tMin = 0;
  _tMax = 1;
  _dt = 1;

  _mode = Linear;

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

inline irtkTemporalHomogeneousTransformation::irtkTemporalHomogeneousTransformation(double dt, double t1, double t2)
{
  int i;

  // Initialize control point parameters (time)
  _tMin = t1;
  _tMax = t2;
  _Nt = round((_tMax - _tMin) / dt) + 1;
  _dt = (_tMax - _tMin) / (_Nt - 1);

  _mode = Linear;

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

inline irtkTemporalHomogeneousTransformation::irtkTemporalHomogeneousTransformation(const irtkMatrix *matrix, double dt, double t1, double t2)
{
  int i;

  // Initialize control point parameters (time)
  _tMin = t1;
  _tMax = t2;
  _Nt = round((_tMax - _tMin) / dt) + 1;
  _dt = (_tMax - _tMin) / (_Nt - 1);

  _mode = Linear;

  _matrix = this->Allocate(_matrix, _Nt);
//  _matrix = new irtkMatrix[_Nt];
  for (i = 0; i < _Nt; i++) {
    _matrix[i] = matrix[i];
  }
//  for (i = -4; i < 0; i++) {
//    _matrix[i] = matrix[0];
//  }
//  for (i = _Nt; i < _Nt+4; i++) {
//    _matrix[i] = matrix[_Nt-1];
//  }

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = _Active;
  }

  _transformMatrix = _matrix[0];
  _tCurrent = _tMin;
}

inline irtkTemporalHomogeneousTransformation::irtkTemporalHomogeneousTransformation(const irtkTemporalHomogeneousTransformation &t) : irtkTransformation(t)
{
  int i;

  _tMin = t._tMin;
  _tMax = t._tMax;
  _Nt = t._Nt;
  _dt = t._dt;

  _mode = Linear;

  _matrix = this->Allocate(_matrix, _Nt);
//  _matrix = new irtkMatrix[_Nt];
  for (i = 0; i < _Nt; i++) {
    _matrix[i] = t._matrix[i];
  }
//  for (i = -4; i < 0; i++) {
//    _matrix[i] = matrix[0];
//  }
//  for (i = _Nt; i < _Nt+4; i++) {
//    _matrix[i] = matrix[_Nt-1];
//  }

  // Allocate memory for DOF status
  _status = new _Status[t.NumberOfDOFs()];
  // Initialize memory for DOF status
  for (i = 0; i < t.NumberOfDOFs(); i++) {
    _status[i] = t._status[i];
  }

  _transformMatrix = _matrix[0];
  _tCurrent = _tMin;
}

inline irtkTemporalHomogeneousTransformation::~irtkTemporalHomogeneousTransformation()
{
  // Free memory for control points if necessary
  if (_matrix != NULL) _matrix = this->Deallocate(_matrix, _Nt);

  _Nt = 0;
}

inline void irtkTemporalHomogeneousTransformation::Initialize(double dt, double t1, double t2)
{
  int i;

  // Free memory for control points if necessary
  if (_matrix != NULL) _matrix = this->Deallocate(_matrix, _Nt);

  // Initialize control point parameters (time)
  _tMin = t1;
  _tMax = t2;
  _Nt = round((_tMax - _tMin) / dt) + 1;
  _dt = (_tMax - _tMin) / (_Nt - 1);

  _mode = Linear;

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

inline void irtkTemporalHomogeneousTransformation::PutMatrix(const irtkMatrix &matrix, int i)
{
  if (i < 0 || i >= _Nt)
	cerr << "wrong index "<<i<<endl;

  _matrix[i] = matrix;
}

inline void irtkTemporalHomogeneousTransformation::PutMatrix(irtkMatrix *matrix)
{
  _matrix = matrix;
}

inline void irtkTemporalHomogeneousTransformation::GlobalTransform(double &x, double &y, double &z, double t)
{
  this->Transform(x, y, z, t);
}

inline void irtkTemporalHomogeneousTransformation::LocalTransform(double &, double &, double &, double)
{
}

inline void irtkTemporalHomogeneousTransformation::GlobalDisplacement(double &x, double &y, double &z, double t)
{
  double a, b, c;

  a = x;
  b = y;
  c = z;
  this->Transform(a, b, c, t);
  x = a - x;
  y = b - y;
  z = c - z;
}

inline void irtkTemporalHomogeneousTransformation::LocalDisplacement(double &x, double &y, double &z, double)
{
  x = 0;
  y = 0;
  z = 0;
}

inline double irtkTemporalHomogeneousTransformation::LatticeToTime(double t) const
{
  return t*(_tMax - _tMin)/double(_Nt - 1)+_tMin;
}

inline double irtkTemporalHomogeneousTransformation::TimeToLattice(double t) const
{
  return (t - _tMin)*(_Nt - 1)/(_tMax - _tMin);
}

inline void irtkTemporalHomogeneousTransformation::SetIntMode(int mode)
{
  if (mode < Linear || mode > CSpline) {
	cerr << "mode " << mode << " not defined" << endl;
	cerr << "use Linear or CSpline" << endl;
	exit(1);
  }
  if (mode == CSpline && _Nt < 4) {
	cerr << "mode CSpline makes no sense with _Nt < 4" << endl;
	exit(1);
  }
  _mode = mode;
}

inline int irtkTemporalHomogeneousTransformation::GetInMode() const
{
  return _mode;
}

inline const char *irtkTemporalHomogeneousTransformation::NameOfClass()
{
  return "irtkTemporalHomogeneousTransformation";
}

#endif
