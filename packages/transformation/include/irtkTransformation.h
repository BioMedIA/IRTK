/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKTRANSFORMATION_H

#define _IRTKTRANSFORMATION_H

#include <irtkImage.h>
#include <irtkGeometry.h>

/// Definition of available states for individual degrees of freedom
typedef enum { _Active, _Passive, _Unknown } _Status;

#define IRTKTRANSFORMATION_MAGIC            815007

#define IRTKTRANSFORMATION_HOMOGENEOUS      1
#define IRTKTRANSFORMATION_RIGID            2
#define IRTKTRANSFORMATION_AFFINE           3
#define IRTKTRANSFORMATION_BSPLINE_FFD      4
#define IRTKTRANSFORMATION_LINEAR_FFD       5
#define IRTKTRANSFORMATION_EIGEN_FFD        6
#define IRTKTRANSFORMATION_MFFD             7
#define IRTKTRANSFORMATION_FLUID            8
#define IRTKTRANSFORMATION_LATTICE_FFD      9
#define IRTKTRANSFORMATION_MULTI_FRAME_LATTICE_FFD 10
#define IRTKTRANSFORMATION_QUATERNION       11
#define IRTKTRANSFORMATION_BSPLINE_FFD_EXT1 12
#define IRTKTRANSFORMATION_LINEAR_FFD_EXT1  13
#define IRTKTRANSFORMATION_BSPLINE_FFD_4D   14
#define IRTKTRANSFORMATION_PERIODIC         20

#define FFDLOOKUPTABLESIZE 1000


/**
 * Abstract base class for general transformations.
 *
 * This is the abstract base class which defines a common interface for all
 * transformation. Each derived class has to implement all abstract member
 * functions.
 */

class irtkTransformation : public irtkObject
{

protected:

  /// Status of each degree of freedom (Active or Passive)
  _Status *_status;

  /// Constructor sets _status to NULL on behalf of subclasses
  irtkTransformation();

public:

  /** Static constructor. This functions returns a pointer to a concrete
   *  transformation by reading the transformation parameters from a file
   *  and creating the approriate transformation
   */
  static irtkTransformation *New(char *);

  /** Static constructor. This functions returns a pointer to a concrete
   *  transformation by copying the transformation passed to it
   */
  static irtkTransformation *New(irtkTransformation *);

  /// Virtual destructor
  virtual ~irtkTransformation();

  /// Returns the number of parameters of the transformation (abstract)
  virtual int    NumberOfDOFs() const = 0;

  /// Gets a transformation parameter (abstract)
  virtual double Get(int) const = 0;

  /// Puts a transformation paramater (abstract)
  virtual void   Put(int, double) = 0;

  /// Puts a control point status
  virtual void   PutStatus(int, _Status);

  /// Gets a control point status
  virtual _Status GetStatus(int) const;

  /// Transforms a single point
  virtual void Transform(irtkPoint &);

  /// Transforms a set of points
  virtual void Transform(irtkPointSet &);

  /// Transforms a single point in 4D
  virtual void Transform(double &, double &, double &, double = 0) = 0;

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0) = 0;

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double &, double &, double &, double = 0) = 0;

  /// Calculate displacement
  virtual void Displacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &, double = 0);

  /// Inverts the transformation (abstract)
  virtual double Inverse(double &, double &, double &, double = 0, double = 0.01) = 0;

  /// Calculate the Jacobian of the transformation with respect to the transformation parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0);

  /// Calculate the Jacobian of the transformation with respect to world coordinates
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0) = 0;

  /// Calculate the Jacobian of the local transformation with respect to world coordinates
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0) = 0;

  /// Calculate the Jacobian of the global transformation with respect to world coordinates
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0) = 0;

  /// Calculate the determinant of the Jacobian of the transformation with respect to world coordinates
  double Jacobian(double, double, double, double = 0);

  /// Calculate the determinant of the Jacobian of the local transformation with respect to world coordinates
  double LocalJacobian(double, double, double, double = 0);

  /// Calculate the determinant of the Jacobian of the global transformation with respect to world coordinates
  double GlobalJacobian(double, double, double, double = 0);

  /// Calculate displacement vectors for image
  virtual void Displacement(irtkGenericImage<double> &, double = 0);

  /// Checks whether transformation is an identity mapping (abstract)
  virtual bool IsIdentity() = 0;

  /// Reads a transformation from a file
  virtual void Read (char *);

  /// Writes a transformation to a file
  virtual void Write(char *);

  /// Imports a transformation from a file
  virtual void Import(char *);

  /// Exports a transformation to a file
  virtual void Export(char *);

  /// Reads a transformation from a file (abstract)
  virtual irtkCifstream& Read(irtkCifstream&) = 0;

  /// Writes a transformation to a file (abstract)
  virtual irtkCofstream& Write( irtkCofstream&) = 0;

  /// Imports a transformation from a file
  virtual istream& Import(istream&);

  /// Exports a transformation to a file
  virtual ostream& Export(ostream&);

  /// I/O
  virtual void Draw ();

  /// Prints the parameters of the transformation (abstract)
  virtual void Print() = 0;

  /// Returns a string with the name of the instantiated class (abstract)
  virtual const char *NameOfClass() = 0;
};


inline _Status irtkTransformation::GetStatus(int index) const
{
  if (index >= this->NumberOfDOFs()) {
    cerr << "irtkTransformation::GetStatus: No such dof" << endl;
    exit(1);
  }
  return _status[index];
}

inline void irtkTransformation::PutStatus(int index, _Status s)
{
  if (index >= this->NumberOfDOFs()) {
    cerr << "irtkTransformation::PutStatus: No such dof" << endl;
    exit(1);
  }
  _status[index] = s;
}

inline void irtkTransformation::Transform(irtkPoint &p)
{
  this->Transform(p._x, p._y, p._z);
}

inline void irtkTransformation::Transform(irtkPointSet &pset)
{
  for (int i = 0; i < pset.Size(); i++) {
    irtkPoint p = pset(i);
    this->Transform(p._x, p._y, p._z);
    pset(i) = p;
  }
}

inline void irtkTransformation::Displacement(double &x, double &y, double &z, double t)
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

inline void irtkTransformation::GlobalDisplacement(double &x, double &y, double &z, double t)
{
  double a, b, c;

  a = x;
  b = y;
  c = z;
  this->GlobalTransform(a, b, c, t);
  x = a - x;
  y = b - y;
  z = c - z;
}

inline void irtkTransformation::LocalDisplacement(double &x, double &y, double &z, double t)
{
  double a, b, c;

  a = x;
  b = y;
  c = z;
  this->LocalTransform(a, b, c, t);
  x = a - x;
  y = b - y;
  z = c - z;
}

inline void irtkTransformation::JacobianDOFs(double [3], int, double, double, double, double)
{
	cerr << this->NameOfClass() << ": JacobianDOFs() not implemented for this class" << endl;
	exit(1);
}

inline double irtkTransformation::Jacobian(double x, double y, double z, double t)
{
  irtkMatrix jac(3, 3);

  // Calculate Jacobian
  this->Jacobian(jac, x, y, z, t);

  // Determinant of Jacobian of deformation derivatives
  return (jac(0, 0)*jac(1, 1)*jac(2, 2) + jac(0, 1)*jac(1, 2)*jac(2, 0) +
          jac(0, 2)*jac(1, 0)*jac(2, 1) - jac(0, 2)*jac(1, 1)*jac(2, 0) -
          jac(0, 0)*jac(1, 2)*jac(2, 1) - jac(0, 1)*jac(1, 0)*jac(2, 2));
}


inline double irtkTransformation::LocalJacobian(double x, double y, double z, double t)
{
  irtkMatrix jac(3, 3);

  // Calculate Jacobian
  this->LocalJacobian(jac, x, y, z, t);

  // Determinant of Jacobian of deformation derivatives
  return (jac(0, 0)*jac(1, 1)*jac(2, 2) + jac(0, 1)*jac(1, 2)*jac(2, 0) +
          jac(0, 2)*jac(1, 0)*jac(2, 1) - jac(0, 2)*jac(1, 1)*jac(2, 0) -
          jac(0, 0)*jac(1, 2)*jac(2, 1) - jac(0, 1)*jac(1, 0)*jac(2, 2));
}

inline double irtkTransformation::GlobalJacobian(double x, double y, double z, double t)
{
  irtkMatrix jac(3, 3);

  // Calculate Jacobian
  this->GlobalJacobian(jac, x, y, z, t);

  // Determinant of Jacobian of deformation derivatives
  return (jac(0, 0)*jac(1, 1)*jac(2, 2) + jac(0, 1)*jac(1, 2)*jac(2, 0) +
          jac(0, 2)*jac(1, 0)*jac(2, 1) - jac(0, 2)*jac(1, 1)*jac(2, 0) -
          jac(0, 0)*jac(1, 2)*jac(2, 1) - jac(0, 1)*jac(1, 0)*jac(2, 2));
}

inline void irtkTransformation::Import(char *name)
{
  ifstream from(name);

  if (!from) {
    cerr << "irtkTransformation::Import: Can't open file " << name << "\n";
    exit(1);
  }

  this->Import(from);
}

inline void irtkTransformation::Export(char *name)
{
  ofstream to(name);

  if (!to) {
    cerr << "irtkTransformation::Export: Can't open file " << name << "\n";
    exit(1);
  }

  this->Export(to);

  // This seems to be required under Mac OS X
  to.close();
}

inline istream& irtkTransformation::Import(istream &)
{
  cerr << "irtkTransformation::Import: Cannot import transformation" << endl;
  exit(1);
}

inline ostream& irtkTransformation::Export(ostream &)
{
  cerr << "irtkTransformation::Export: Cannot export transformation" << endl;
  exit(1);
}


// Homogeneous transformation classes
#include <irtkHomogeneousTransformation.h>
#include <irtkRigidTransformation.h>
#include <irtkAffineTransformation.h>

// Free-form transformations
#include <irtkFreeFormTransformation.h>
#include <irtkFreeFormTransformation3D.h>
#include <irtkFreeFormTransformation4D.h>

#include <irtkLinearFreeFormTransformation.h>
#include <irtkBSplineFreeFormTransformation3D.h>
#include <irtkEigenFreeFormTransformation.h>

#include <irtkBSplineFreeFormTransformation4D.h>
#include <irtkBSplineFreeFormTransformationPeriodic.h>

// Composite transformations
#include <irtkMultiLevelFreeFormTransformation.h>
#include <irtkFluidFreeFormTransformation.h>

// Image transformation filters
#include <irtkImageTransformation.h>
#include <irtkImageTransformation2.h>
#include <irtkImageHomogeneousTransformation.h>

// Typedefs for backwards compatibility
typedef class irtkBSplineFreeFormTransformation3D irtkBSplineFreeFormTransformation;


#endif
