





#ifndef _IRTKBSPLINEFREEFORMTRANSFORMATIONPERIODIC_H

#define _IRTKBSPLINEFREEFORMTRANSFORMATIONPERIODIC_H

#include <irtkGeometry.h>

/**
 * Class for free form transformations based on tensor product B-splines.
 *
 * This class implements 4D free form transformation using B-splines.
 *
 * For more details about the implementation see Lee, Wolberg and Shin, IEEE
 * Transactions on Visualization and Computer Graphics, Vol. 3, No. 3, 1997.
 */

class irtkBSplineFreeFormTransformationPeriodic : public irtkFreeFormTransformation4D
{

protected:

  int _periodic;

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

  /// Calculate the bending energy of the transformation at control points (2D)
  virtual double Bending2D(int i, int j, int t);

  /// Calculate the bending energy of the transformation at control points (3D)
  virtual double Bending3D(int i, int j, int k, int t);

  /// Calculate the gradient of the bending energy with respect to the parameters
  virtual void BendingGradient2D(double *gradient);

  /// Calculate the gradient of the bending energy with respect to the parameters
  virtual void BendingGradient3D(double *gradient);

  /// Memory for lookup table for B-spline basis function values
  static    double LookupTable[FFDLOOKUPTABLESIZE][4];

  /* Memory for lookup table for first derivatives of B-spline basis
   * function values */
  static    double LookupTable_I[FFDLOOKUPTABLESIZE][4];

  /* Memory for lookup table for second derivatives of B-spline basis
   * function values */
  static    double LookupTable_II[FFDLOOKUPTABLESIZE][4];

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
  irtkBSplineFreeFormTransformationPeriodic();

  /// Constructor
  irtkBSplineFreeFormTransformationPeriodic(irtkBaseImage &, double, double, double, double);

  /// Constructor
  irtkBSplineFreeFormTransformationPeriodic(double x1, double y1, double z1, double t1,
											double x2, double y2, double z2, double t2,
											double dx, double dy, double dz, double dt,
											double* xaxis, double* yaxis, double* zaxis);

  /// Copy Constructor
  irtkBSplineFreeFormTransformationPeriodic(const class irtkBSplineFreeFormTransformationPeriodic &);

  /// Destructor
  virtual ~irtkBSplineFreeFormTransformationPeriodic();

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *, double *,
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

  /// Subdivide FFD in 2D
  virtual void Subdivide2D();

  /// Subdivide FFD in 3D
  virtual void Subdivide3D();

  /// Subdivide FFD in 4D
  virtual void Subdivide4D();

  /// Calculates the FFD (for a point in FFD coordinates) with checks
  virtual void FFD1(double &, double &, double &, double) const;

  /// Calculates the FFD (for a point in FFD coordinates) without checks
  virtual void FFD2(double &, double &, double &, double) const;

  /// Transforms a point
  virtual void Transform(double &, double &, double &, double = 0);

  /// Transforms a point
  virtual void Transform2(double &, double &, double &, double = 0);

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

  /// new (periodic) TimeToLattice function
  double TimeToLattice(double) const;

  /// Calculate the Jacobian of the transformation
  virtual void Jacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the local transformation
  virtual void LocalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate the Jacobian of the global transformation
  virtual void GlobalJacobian(irtkMatrix &, double, double, double, double = 0);

  /// Calculate total bending energy
  virtual double Bending();

  /// Calculate the bending energy of the transformation
  virtual double Bending(double, double, double, double);

  /// Calculate the gradient of the bending energy with respect to the parameters
  virtual void BendingGradient(double *gradient);

  /** Returns the bounding box for a control point (in mm). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBoxCP(int, irtkPoint &, irtkPoint &, double = 1) const;

  /** Returns the bounding box for a control point. The last parameter
   *  specifies what fraction of the bounding box to return. The default
   *  is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBoxCP(int, double &, double &, double &, double &, double &,
                           double &, double &, double &, double = 1) const;

  /** Returns the bounding box for a control point (in pixels). The last
   *  parameter specifies what fraction of the bounding box to return. The
   *  default is 1 which equals 100% of the bounding box.
   */
  virtual void BoundingBoxImage(irtkGreyImage *, int, int &, int &, int &, int &,
                           int &, int &, int &, int &, double = 1) const;
  // temporal bounding box is in continuous space (no discrete frames!)
  virtual void BoundingBoxImage(irtkGreyImage *, int, int &, int &, int &, int &,
                           int &, int &, double &, double &, double = 1) const;
  // for inheritance
  virtual void BoundingBoxImage(irtkGreyImage *, int, int &, int &, int &, int &,
                           int &, int &, double = 1) const;

//  // alternative bounding box without index
//  virtual void BoundingBox(irtkGreyImage *, int &, int &, int &, int &,
//                           int &, int &, double &, double &, double = 1) const;

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

  /// Use periodic time on/off
  void PeriodicOn();
  void PeriodicOff();

  /// Calculates the 3D FFD at time t
  irtkBSplineFreeFormTransformation3D *Compute3DFFD1(double);
  irtkBSplineFreeFormTransformation3D *Compute3DFFD2(double);
  irtkBSplineFreeFormTransformation3D *Compute3DFFD3(double);

};

inline void irtkBSplineFreeFormTransformationPeriodic::PeriodicOn()
{
  _periodic = true;
}

inline void irtkBSplineFreeFormTransformationPeriodic::PeriodicOff()
{
  _periodic = false;
}

inline double irtkBSplineFreeFormTransformationPeriodic::TimeToLattice(double t) const
{
  double tt = t;
  if (_periodic) {
    while (tt < 0)
      tt += 1.;
    while (tt >= 1)
      tt -= 1.;
  }
  tt = (tt - _tMin)*(_t - 1)/(_tMax - _tMin);
  return tt;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B(double x)
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

inline double irtkBSplineFreeFormTransformationPeriodic::B(int i, double t)
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

inline double irtkBSplineFreeFormTransformationPeriodic::B_I(int i, double t)
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

inline double irtkBSplineFreeFormTransformationPeriodic::B_II(int i, double t)
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

inline double irtkBSplineFreeFormTransformationPeriodic::B0(double t)
{
  return (1-t)*(1-t)*(1-t)/6.0;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B1(double t)
{
  return (3*t*t*t - 6*t*t + 4)/6.0;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B2(double t)
{
  return (-3*t*t*t + 3*t*t + 3*t + 1)/6.0;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B3(double t)
{
  return (t*t*t)/6.0;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B0_I(double t)
{
  return -(1-t)*(1-t)/2.0;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B1_I(double t)
{
  return (9*t*t - 12*t)/6.0;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B2_I(double t)
{
  return (-9*t*t + 6*t + 3)/6.0;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B3_I(double t)
{
  return (t*t)/2.0;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B0_II(double t)
{
  return 1 - t;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B1_II(double t)
{
  return 3*t - 2;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B2_II(double t)
{
  return -3*t + 1;
}

inline double irtkBSplineFreeFormTransformationPeriodic::B3_II(double t)
{
  return t;
}

inline void irtkBSplineFreeFormTransformationPeriodic::Transform(double &x, double &y, double &z, double t)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  this->FFD1(u, v, w, this->TimeToLattice(t));

  // Add FFD to world coordinates
  x += u;
  y += v;
  z += w;
}

inline void irtkBSplineFreeFormTransformationPeriodic::Transform2(double &x, double &y, double &z, double t)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  this->FFD2(u, v, w, this->TimeToLattice(t));

  // Add FFD to world coordinates
  x += u;
  y += v;
  z += w;
}

inline void irtkBSplineFreeFormTransformationPeriodic::GlobalTransform(double &, double &, double &, double)
{}

inline void irtkBSplineFreeFormTransformationPeriodic::LocalTransform(double &x, double &y, double &z, double t)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // calculate displacement
  this->LocalDisplacement(u, v, w, t);

  // Add displacement
  x += u;
  y += v;
  z += w;
}

inline void irtkBSplineFreeFormTransformationPeriodic::GlobalDisplacement(double &x, double &y, double &z, double)
{
  x = 0;
  y = 0;
  z = 0;
}

inline void irtkBSplineFreeFormTransformationPeriodic::LocalDisplacement(double &x, double &y, double &z, double t)
{
  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(x, y, z);

  // Calculate FFD
  this->FFD1(x, y, z, this->TimeToLattice(t));
}

inline void irtkBSplineFreeFormTransformationPeriodic::Displacement(double &x, double &y, double &z, double t)
{
	this->LocalDisplacement(x, y, z, t);
}

inline const char *irtkBSplineFreeFormTransformationPeriodic::NameOfClass()
{
  return "irtkBSplineFreeFormTransformationPeriodic";
}

#endif
