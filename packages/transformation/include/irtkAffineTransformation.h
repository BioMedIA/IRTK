/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKAFFINETRANSFORMATION_H

#define _IRTKAFFINETRANSFORMATION_H

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

class irtkAffineTransformation : public irtkRigidTransformation
{

protected:

  /// Scaling along the x-axis
  double _sx;

  /// Scaling along the y-axis
  double _sy;

  /// Scaling along the z-axis
  double _sz;

  /// Skew angle in the x direction based on y component (in degrees)
  double _sxy;

  /// Skew angle in the y direction based on z component (in degrees)
  double _syz;

  /// Skew angle in the x direction based on z component (in degrees)
  double _sxz;

  /// Construct a matrix based on parameters passed in the array.
  virtual irtkMatrix Parameters2Matrix(double *) const;

  /// Construct a matrix based on old style parameters.
  virtual irtkMatrix Parameters2Matrix_15DOFs(double *) const;

  /// Assign the parameters passed in the array to the current object and
  /// update the matrix.
  virtual void SetParameters(double *);

public:

  /// Constructor (default)
  irtkAffineTransformation();

  /// Reset Transformation
  virtual void Reset();

  /// Constructor (copy)
  irtkAffineTransformation(const irtkRigidTransformation &);

  /// Constructor (copy)
  irtkAffineTransformation(const irtkAffineTransformation &);

  /// Destructor
  virtual ~irtkAffineTransformation();

  /// Puts scaling factor along the x-axis
  virtual void   PutScaleX(double);

  /// Gets scaling factor along the x-axis
  virtual double GetScaleX();

  /// Puts scaling factor along the y-axis
  virtual void   PutScaleY(double);

  /// Gets scaling factor along the y-axis
  virtual double GetScaleY();

  /// Puts scaling factor along the z-axis
  virtual void   PutScaleZ(double);

  /// Gets scaling factor along the z-axis
  virtual double GetScaleZ();

  /// Puts y-dependent skewing angle in the x direction (in degrees)
  virtual void   PutShearXY(double);

  /// Gets y-dependent skewing angle in the x direction (in degrees)
  virtual double GetShearXY();

  /// Puts z-dependent skewing angle in the y direction (in degrees)
  virtual void   PutShearYZ(double);

  /// Gets z-dependent skewing angle in the y direction (in degrees)
  virtual double GetShearYZ();

  /// Puts z-dependent skewing angle in the x direction (in degrees)
  virtual void   PutShearXZ(double);

  /// Gets z-dependent skewing angle in the x direction (in degrees)
  virtual double GetShearXZ();

  /// Returns the number of parameters of the transformation
  virtual int NumberOfDOFs() const;

  /// Puts a transformation parameter (transformation matrix is updated)
  virtual void   Put(int, double);

  /// Gets a transformation parameter
  virtual double Get(int) const;

  /// Checks whether transformation is an identity mapping
  virtual Bool IsIdentity();

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

  /// Puts the transformation matrix (transformation parameters are updated)
  virtual void PutMatrix(const irtkMatrix &);

  /// Update transformation matrix
  virtual void UpdateMatrix();

  /// Updates transformation parameters based on current matrix.
  virtual void UpdateParameter();

};

inline void irtkAffineTransformation::Reset()
{
  // Initialize rotations and translations
  _tx = _ty = _tz = 0;
  _rx = _ry = _rz = 0;

  // Initialize scale and shears
  _sx  = _sy  = _sz  = 100;
  _sxy = _syz = _sxz = 0;

  //define transformation's orientation
  /*irtkMatrix m(3,3); double params[3];
  int RX = 0; int RY = 1; int RZ = 2;
  double *x,*y,*z;
  x = new double[3];
  y = new double[3];
  z = new double[3];
  input.GetOrientation(x,y,z); 
  m(0,0) = x[0]; m(0,1) = x[1]; m(0,2) = x[2];
  m(1,0) = y[0]; m(1,1) = y[1]; m(1,2) = y[2];
  m(2,0) = z[0]; m(2,1) = z[1]; m(2,2) = z[2];  
  delete x,y,z;
  double tmp = asin(-1 * m(0,2));
  //double tmp1 = m(0,2);
  //double cost = cos(tmp);

  // asin returns values for tmp in range -pi/2 to +pi/2, i.e. cos(tmp) >=
  // 0 so the division by cos(tmp) in the first part of the if clause was
  // not needed.
  if (fabs(cos(tmp)) > 0.0000001) {
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

  //set rotation for transformation
  _rx = params[RX];
  _ry = params[RY];
  _rz = params[RZ];*/
  // Update transformation matrix
  this->UpdateMatrix();
}


inline int irtkAffineTransformation::NumberOfDOFs() const
{
  return 12;
}

inline irtkAffineTransformation::irtkAffineTransformation()
{
  int i;

  // Initialize rotations and translations
  _tx = _ty = _tz = 0;
  _rx = _ry = _rz = 0;

  // Initialize scale and shears
  _sx  = _sy  = _sz  = 100;
  _sxy = _syz = _sxz = 0;

  // Free memory allocated for DOF status by base class
  delete []_status;

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = _Active;
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

inline irtkAffineTransformation::irtkAffineTransformation(const irtkRigidTransformation &t)
{
  int i;

  // Copy rotations and translations
  _tx = t.Get(0);
  _ty = t.Get(1);
  _tz = t.Get(2);
  _rx = t.Get(3);
  _ry = t.Get(4);
  _rz = t.Get(5);

  // Initialize scale and shears
  _sx  = 100;
  _sy  = 100;
  _sz  = 100;
  _sxy = 0;
  _syz = 0;
  _sxz = 0;

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

  // Update transformation matrix
  this->UpdateMatrix();
}

inline irtkAffineTransformation::irtkAffineTransformation(const irtkAffineTransformation &t)
{
  int i;

  // Copy rotations and translations
  _tx = t._tx;
  _ty = t._ty;
  _tz = t._tz;
  _rx = t._rx;
  _ry = t._ry;
  _rz = t._rz;

  // Copy scale and shears
  _sx  = t._sx;
  _sy  = t._sy;
  _sz  = t._sz;
  _sxy = t._sxy;
  _syz = t._syz;
  _sxz = t._sxz;

  // Free memory allocated for DOF status by base class
  delete []_status;

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = t._status[i];
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

inline irtkAffineTransformation::~irtkAffineTransformation()
{

}

inline void irtkAffineTransformation::PutScaleX(double sx)
{
  _sx = sx;
  this->UpdateMatrix();
}

inline double irtkAffineTransformation::GetScaleX()
{
  return _sx;
}

inline void irtkAffineTransformation::PutScaleY(double sy)
{
  _sy = sy;
  this->UpdateMatrix();
}

inline double irtkAffineTransformation::GetScaleY()
{
  return _sy;
}

inline void irtkAffineTransformation::PutScaleZ(double sz)
{
  _sz = sz;
  this->UpdateMatrix();
}

inline double irtkAffineTransformation::GetScaleZ()
{
  return _sz;
}

inline void irtkAffineTransformation::PutShearXY(double sxy)
{
  _sxy = sxy;
  this->UpdateMatrix();
}

inline double irtkAffineTransformation::GetShearXY()
{
  return _sxy;
}

inline void irtkAffineTransformation::PutShearYZ(double syz)
{
  _syz = syz;
  this->UpdateMatrix();
}

inline double irtkAffineTransformation::GetShearYZ()
{
  return _syz;
}

inline void irtkAffineTransformation::PutShearXZ(double sxz)
{
  _sxz = sxz;
  this->UpdateMatrix();
}

inline double irtkAffineTransformation::GetShearXZ()
{
  return _sxz;
}

inline const char *irtkAffineTransformation::NameOfClass()
{
  return "irtkAffineTransformation";
}

inline void irtkAffineTransformation::PutMatrix(const irtkMatrix &matrix)
{
  _matrix = matrix;
  this->UpdateParameter();
}

inline void irtkAffineTransformation::SetParameters(double *params)
{
  this->irtkRigidTransformation::SetParameters(params);

  _sx  = params[SX];
  _sy  = params[SY];
  _sz  = params[SZ];

  _sxy = params[SXY];
  _syz = params[SYZ];
  _sxz = params[SXZ];

  this->UpdateMatrix();
}

#endif
