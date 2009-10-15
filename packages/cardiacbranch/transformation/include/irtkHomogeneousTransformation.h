/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _HOMOGENEOUSTRANSFORMATION_H

#define _HOMOGENEOUSTRANSFORMATION_H

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

typedef enum {TX, TY, TZ, RX, RY, RZ, SX, SY, SZ, SXY, SYZ, SXZ, SYX, SZY, SZX}
irtkHomogeneousTransformationParameterIndex;

class irtkHomogeneousTransformation : public irtkTransformation
{

protected:

  /// 4 x 4 transformation matrix for homogeneous coordinates
  irtkMatrix _matrix;

public:

  /// Constructor (default)
  irtkHomogeneousTransformation();

  /// Constructor (from matrix)
  irtkHomogeneousTransformation(const irtkMatrix &);

  /// Constructor (copy)
  irtkHomogeneousTransformation(const irtkHomogeneousTransformation &);

  /// Destructor
  virtual ~irtkHomogeneousTransformation();

  /// Returns the number of parameters of the transformation
  virtual int    NumberOfDOFs() const;

  /// Gets a transformation parameter
  virtual double Get(int) const;

  /// Puts a transformation paramater
  virtual void   Put(int, double);

  /// Gets the transformation matrix
  virtual irtkMatrix GetMatrix() const;

  /// Puts the transformation matrix (Argument is not checked)
  virtual void      PutMatrix(const irtkMatrix &);

  /// Transforms a single point
  virtual void Transform(double &, double &, double &, double = 0);

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0);

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double &, double &, double &, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &, double = 0);

  /// Inverts the transformation
  virtual void Invert();

  /// Inverse transformation
  virtual double Inverse(double &, double &, double &, double = 0, double = 0.01);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the global transformation
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Checks whether transformation is an identity mapping
  virtual Bool IsIdentity();

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

inline int irtkHomogeneousTransformation::NumberOfDOFs() const
{
  return 12;
}

inline irtkHomogeneousTransformation::irtkHomogeneousTransformation()
{
  int i;

  // Create 4 x 4 transformation matrix
  _matrix = irtkMatrix(4, 4);

  // Initialize to identity
  _matrix.Ident();

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = _Active;
  }
}

inline irtkHomogeneousTransformation::irtkHomogeneousTransformation(const irtkMatrix &matrix)
{
  int i;

  _matrix = matrix;

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = _Active;
  }
}

inline irtkHomogeneousTransformation::irtkHomogeneousTransformation(const irtkHomogeneousTransformation &t)
{
  int i;

  _matrix = t._matrix;

  // Allocate memory for DOF status
  _status = new _Status[this->NumberOfDOFs()];

  // Initialize memory for DOF status
  for (i = 0; i < this->NumberOfDOFs(); i++) {
    _status[i] = t._status[i];
  }
}

inline irtkHomogeneousTransformation::~irtkHomogeneousTransformation()
{
}

inline void irtkHomogeneousTransformation::PutMatrix(const irtkMatrix &matrix)
{
  _matrix = matrix;
}

inline irtkMatrix irtkHomogeneousTransformation::GetMatrix() const
{
  return _matrix;
}

inline void irtkHomogeneousTransformation::GlobalTransform(double &x, double &y, double &z, double t)
{
  this->Transform(x, y, z, t);
}

inline void irtkHomogeneousTransformation::LocalTransform(double &x, double &y, double &z, double)
{
}

inline void irtkHomogeneousTransformation::GlobalDisplacement(double &x, double &y, double &z, double t)
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

inline void irtkHomogeneousTransformation::LocalDisplacement(double &x, double &y, double &z, double)
{
  x = 0;
  y = 0;
  z = 0;
}

inline const char *irtkHomogeneousTransformation::NameOfClass()
{
  return "irtkHomogeneousTransformation";
}

#endif
