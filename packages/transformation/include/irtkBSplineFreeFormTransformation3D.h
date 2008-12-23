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

  /** Computes the B-spline coefficients need to interpolate a displacement
      field. */
  void ComputeCoefficients(double* dxs, double* dys, double* dzs, irtkRealImage& xCoeffs, irtkRealImage& yCoeffs, irtkRealImage& zCoeffs);

  /// Memory for lookup table for B-spline basis function values
  static    double LookupTable[FFDLOOKUPTABLESIZE][4];

  /* Memory for lookup table for first derivatives of B-spline basis
   * function values */
  static    double LookupTable_I[FFDLOOKUPTABLESIZE][4];

  /* Memory for lookup table for second derivatives of B-spline basis
   * function values */
  static    double LookupTable_II[FFDLOOKUPTABLESIZE][4];

public:

  /// Returns the value of the i-th B-spline basis function
  static double B (int, double);

  /// Returns the 1st derivative value of the i-th B-spline basis function
  static double B_I(int, double);

  /// Returns the 2nd derivative value of the i-th B-spline basis function
  static double B_II(int, double);

  /// Constructor
  irtkBSplineFreeFormTransformation3D();

  /// Constructor
  irtkBSplineFreeFormTransformation3D(irtkBaseImage &, double = 1, double = 1, double = 1);

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
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *,
                             double *, double *, double *, int);

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dxs The x-displacements at each control point.
      \param dys The y-displacements at each control point.
      \param dzs The z-displacements at each control point. */
  virtual void Interpolate(double* dxs, double* dys, double* dzs);

  /// Subdivide FFD
  virtual void Subdivide();

  /// Calculates the FFD (for a point in FFD coordinates) with checks
  virtual void FFD1(double &, double &, double &) const;

  /// Calculates the FFD (for a point in FFD coordinates) without checks
  virtual void FFD2(double &, double &, double &) const;

  /// Transforms a point
  virtual void Transform(double &, double &, double &, double = 0);

  /// Transforms a point without checks
  virtual void Transform2(double &, double &, double &, double = 0);

  /// Transforms a point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0);

  /// Transforms a point using the local transformation component only
  virtual void LocalTransform (double &, double &, double &, double = 0);

  /// Calculates displacement using the global transformation component only
  virtual void GlobalDisplacement(double &, double &, double &, double = 0);

  /// Calculates displacement using the local transformation component only
  virtual void LocalDisplacement(double &, double &, double &, double = 0);

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the global transformation
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the bending energy of the transformation
  virtual double Bending(double x, double y, double z);

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
  this->FFD1(u, v, w);

  // Add FFD to world coordinates
  x += u;
  y += v;
  z += w;
}

inline void irtkBSplineFreeFormTransformation3D::Transform2(double &x, double &y, double &z, double)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  this->FFD2(u, v, w);

  // Add FFD to world coordinates
  x += u;
  y += v;
  z += w;
}

inline void irtkBSplineFreeFormTransformation3D::GlobalTransform(double &x, double &y, double &z, double)
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
  this->FFD1(x, y, z);
}

inline const char *irtkBSplineFreeFormTransformation3D::NameOfClass()
{
  return "irtkBSplineFreeFormTransformation3D";
}

#endif
