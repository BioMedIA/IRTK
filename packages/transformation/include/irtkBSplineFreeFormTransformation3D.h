/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKBSPLINEFREEFORMTRANSFORMATION3D_H

#define _IRTKBSPLINEFREEFORMTRANSFORMATION3D_H

#include <irtkGeometry.h>

/**
 * Class for free form transformations based on tensor product B-splines.
 *
 * This class implements 3D free form transformation using B-splines.
 *
 * For more details about the implementation see Lee, Wolberg and Shin, IEEE
 * Transactions on Visualization and Computer Graphics, Vol. 3, No. 3, 1997.
 */

class irtkBSplineFreeFormTransformation3D : public irtkFreeFormTransformation3D
{

  friend class irtkLinearFreeFormTransformation;

protected:

  /// Returns the value of the first B-spline basis function
  static double B0(double);

  /// Returns the value of the second B-spline basis function
  static double B1(double);

  /// Returns the value of the third B-spline basis function
  static double B2(double);

  /// Returns the value of the fourth B-spline basis function
  static double B3(double);

  /// Returns the 1st derivative value of the first B-spline basis function
  static double B0_I(double);

  /// Returns the 1st derivative value of the second B-spline basis function
  static double B1_I(double);

  /// Returns the 1st derivative value of the third B-spline basis function
  static double B2_I(double);

  /// Returns the 1st derivative value of the fourth B-spline basis function
  static double B3_I(double);

  /// Returns the 2nd derivative value of the first B-spline basis function
  static double B0_II(double);

  /// Returns the 2nd derivative value of the second B-spline basis function
  static double B1_II(double);

  /// Returns the 2nd derivative value of the third B-spline basis function
  static double B2_II(double);

  /// Returns the 2nd derivative value of the fourth B-spline basis function
  static double B3_II(double);

  /// Subdivide FFD in 2D
  virtual void Subdivide2D();

  /// Subdivide FFD in 3D
  virtual void Subdivide3D();

  /// Initialize anti-causal coefficients
  static double InitialAntiCausalCoefficient(double *, int, double z);

  /// Initialize causal coefficients
  static double InitialCausalCoefficient(double *, int, double, double);

  /** Convert from displacement field values to B-spline coefficient values. */
  static void ConvertToInterpolationCoefficients(double* c, int DataLength, double* z, int NbPoles, double Tolerance);

  /** Computes the B-spline coefficients needed to interpolate a 2D displacement
      field. */
  void ComputeCoefficients2D(const double* dxs, const double* dys, const double* dzs, irtkGenericImage<double>& xCoeffs, irtkGenericImage<double>& yCoeffs, irtkGenericImage<double>& zCoeffs);

  /** Computes the B-spline coefficients needed to interpolate a 3D displacement
      field. */
  void ComputeCoefficients3D(const double* dxs, const double* dys, const double* dzs, irtkGenericImage<double>& xCoeffs, irtkGenericImage<double>& yCoeffs, irtkGenericImage<double>& zCoeffs);

  /** Approximate displacements in 2D: This function takes a set of points
      and a set of displacements and find a FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate2D(const double *, const double *, const double *, double *, double *, double *, int);

  /** Approximate displacements in 3D: This function takes a set of points
      and a set of displacements and find a FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate3D(const double *, const double *, const double *, double *, double *, double *, int);

  /** Approximate displacements in 2D: This function takes a set of points
      and a set of displacements and find a new FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual void ApproximateAsNew2D(const double *, const double *, const double *, double *, double *, double *, int);

  /** Approximate displacements in 3D: This function takes a set of points
      and a set of displacements and find a new FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual void ApproximateAsNew3D(const double *, const double *, const double *, double *, double *, double *, int);

  /// Memory for lookup table for B-spline basis function values
  static    double LookupTable[FFDLOOKUPTABLESIZE][4];

  /* Memory for lookup table for first derivatives of B-spline basis
   * function values */
  static    double LookupTable_I[FFDLOOKUPTABLESIZE][4];

  /* Memory for lookup table for second derivatives of B-spline basis
   * function values */
  static    double LookupTable_II[FFDLOOKUPTABLESIZE][4];

  /// Calculates a 2D FFD (for a point in FFD coordinates)
  virtual void FFD2D(double &, double &) const;

  /// Calculates a 3D FFD (for a point in FFD coordinates)
  virtual void FFD3D(double &, double &, double &) const;

  /// Calculate the bending energy of the transformation at control points (2D)
  virtual double Bending2D(int i, int j);

  /// Calculate the bending energy of the transformation at arbitrary points  (2D)
  virtual double Bending2D(double x, double y);

  /// Calculate the bending energy of the transformation at control points (3D)
  virtual double Bending3D(int i, int j, int k);

  /// Calculate the bending energy of the transformation at arbitrary points  (3D)
  virtual double Bending3D(double x, double y, double z);

  /// Calculate the gradient of the bending energy with respect to the parameters
  virtual void BendingGradient2D(double *gradient);

  /// Calculate the gradient of the bending energy with respect to the parameters
  virtual void BendingGradient3D(double *gradient);

public:

  /// Returns the value of the B-spline basis function
  static double B (double);

  /// Returns the value of the i-th B-spline basis function
  static double B (int, double);

  /// Returns the 1st derivative value of the i-th B-spline basis function
  static double B_I(int, double);

  /// Returns the 2nd derivative value of the i-th B-spline basis function
  static double B_II(int, double);

  /// Constructor
  irtkBSplineFreeFormTransformation3D();

  /// Constructor
  irtkBSplineFreeFormTransformation3D(irtkBaseImage &, double, double, double);

  /// Constructor
  irtkBSplineFreeFormTransformation3D(double x1, double y1, double z1,
                                      double x2, double y2, double z2,
                                      double dx, double dy, double dz,
                                      double* xaxis, double* yaxis, double* zaxis);

  /// Copy Constructor
  irtkBSplineFreeFormTransformation3D(const class irtkBSplineFreeFormTransformation3D &);

  /// Destructor
  virtual ~irtkBSplineFreeFormTransformation3D();

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(const double *, const double *, const double *, double *, double *, double *, int);

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a !new! FFD which approximates these
      displacements. After approximation the displacements replaced by
      the residual displacement errors at the points */
  virtual void ApproximateAsNew(const double *, const double *, const double *, double *, double *, double *, int);

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dxs The x-displacements at each control point.
      \param dys The y-displacements at each control point.
      \param dzs The z-displacements at each control point. */
  virtual void Interpolate(const double* dxs, const double* dys, const double* dzs);

  /// Subdivide FFD
  virtual void Subdivide();

  /// Transforms a point
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

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void JacobianDetDerivative(irtkMatrix *, int, int, int);

  /// Calculate the Jacobian of the global transformation
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the transformation with respect to the transformation parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0);

  /// Calculate total bending energy
  virtual double Bending();

  /// Calculate bending energy
  virtual double Bending(double, double, double, double = 0);

  /// Calculate the gradient of the bending energy with respect to the parameters
  virtual void BendingGradient(double *gradient);

  /** Returns the bounding box for a control point (in mm). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(int, irtkPoint &, irtkPoint &, double = 1) const;

  /** Returns the bounding box for a control point (in mm). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(int, double &, double &, double &,
                           double &, double &, double &, double = 1) const;

  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBox(irtkGreyImage *, int, int &, int &, int &,
                           int &, int &, int &, double = 1) const;

  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void MultiBoundingBox(irtkGreyImage *, int, int &, int &, int &,
                                int &, int &, int &, double = 1) const;

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

inline double irtkBSplineFreeFormTransformation3D::B(double x)
{
  x = fabs(x);
  double value=0.0;
  if (x < 2.0) {
    if (x < 1.0) {
      value = (double)(2.0f/3.0f + (0.5f*x-1.0)*x*x);
    } else {
      x -=2.0f;
      value = -x*x*x/6.0f;
    }
  }
  return value;
}

inline double irtkBSplineFreeFormTransformation3D::B(int i, double t)
{
  switch (i) {
    case 0:
      return (1-t)*(1-t)*(1-t)/6.0;
    case 1:
      return (3*t*t*t - 6*t*t + 4)/6.0;
    case 2:
      return (-3*t*t*t + 3*t*t + 3*t + 1)/6.0;
    case 3:
      return (t*t*t)/6.0;
  }
  return 0;
}

inline double irtkBSplineFreeFormTransformation3D::B_I(int i, double t)
{
  switch (i) {
    case 0:
      return -(1-t)*(1-t)/2.0;
    case 1:
      return (9*t*t - 12*t)/6.0;
    case 2:
      return (-9*t*t + 6*t + 3)/6.0;
    case 3:
      return (t*t)/2.0;
  }
  return 0;
}

inline double irtkBSplineFreeFormTransformation3D::B_II(int i, double t)
{
  switch (i) {
    case 0:
      return 1 - t;
    case 1:
      return 3*t - 2;
    case 2:
      return -3*t + 1;
    case 3:
      return t;
  }
  return 0;
}

inline double irtkBSplineFreeFormTransformation3D::B0(double t)
{
  return (1-t)*(1-t)*(1-t)/6.0;
}

inline double irtkBSplineFreeFormTransformation3D::B1(double t)
{
  return (3*t*t*t - 6*t*t + 4)/6.0;
}

inline double irtkBSplineFreeFormTransformation3D::B2(double t)
{
  return (-3*t*t*t + 3*t*t + 3*t + 1)/6.0;
}

inline double irtkBSplineFreeFormTransformation3D::B3(double t)
{
  return (t*t*t)/6.0;
}

inline double irtkBSplineFreeFormTransformation3D::B0_I(double t)
{
  return -(1-t)*(1-t)/2.0;
}

inline double irtkBSplineFreeFormTransformation3D::B1_I(double t)
{
  return (9*t*t - 12*t)/6.0;
}

inline double irtkBSplineFreeFormTransformation3D::B2_I(double t)
{
  return (-9*t*t + 6*t + 3)/6.0;
}

inline double irtkBSplineFreeFormTransformation3D::B3_I(double t)
{
  return (t*t)/2.0;
}

inline double irtkBSplineFreeFormTransformation3D::B0_II(double t)
{
  return 1 - t;
}

inline double irtkBSplineFreeFormTransformation3D::B1_II(double t)
{
  return 3*t - 2;
}

inline double irtkBSplineFreeFormTransformation3D::B2_II(double t)
{
  return -3*t + 1;
}

inline double irtkBSplineFreeFormTransformation3D::B3_II(double t)
{
  return t;
}

inline void irtkBSplineFreeFormTransformation3D::Transform(double &x, double &y, double &z, double)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  if (_z == 1) {
  	this->FFD2D(u, v);
  	w = 0;
  } else {
  	this->FFD3D(u, v, w);
  }

  // Add FFD to world coordinates
  x += u;
  y += v;
  z += w;
}

inline void irtkBSplineFreeFormTransformation3D::GlobalTransform(double &, double &, double &, double)
{}

inline void irtkBSplineFreeFormTransformation3D::LocalTransform(double &x, double &y, double &z, double)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // calculate displacement
  this->LocalDisplacement(u, v, w);

  // Add displacement
  x += u;
  y += v;
  z += w;
}

inline void irtkBSplineFreeFormTransformation3D::GlobalDisplacement(double &x, double &y, double &z, double)
{
  x = 0;
  y = 0;
  z = 0;
}

inline void irtkBSplineFreeFormTransformation3D::LocalDisplacement(double &x, double &y, double &z, double)
{

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(x, y, z);

  // Calculate FFD
  if (_z == 1) {
  	this->FFD2D(x, y);
  } else {
  	this->FFD3D(x, y, z);
  }
}

inline void irtkBSplineFreeFormTransformation3D::Displacement(double &x, double &y, double &z, double)
{
  this->LocalDisplacement(x, y, z);
}

inline const char *irtkBSplineFreeFormTransformation3D::NameOfClass()
{
  return "irtkBSplineFreeFormTransformation3D";
}

inline double irtkBSplineFreeFormTransformation3D::Bending(double x, double y, double z, double)
{
  return this->Bending3D(x, y, z);
}

#endif
