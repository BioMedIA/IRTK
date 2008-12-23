/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKBSPLINEBIASFIELD_H

#define _IRTKBSPLINEBIASFIELD_H

#include <irtkGeometry.h>

#define BIASLOOKUPTABLESIZE 1000

class irtkBSplineBiasField : public irtkBiasField
{

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

  /** Computes the B-spline coefficients needed to interpolate a bias
      field. */
  void ComputeCoefficients(double* dbias, irtkRealImage& coeffs);

  /// Memory for lookup table for B-spline basis function values
  static    double LookupTable[BIASLOOKUPTABLESIZE][4];

  /* Memory for lookup table for first derivatives of B-spline basis
   * function values */
  static    double LookupTable_I[BIASLOOKUPTABLESIZE][4];

  /* Memory for lookup table for second derivatives of B-spline basis
   * function values */
  static    double LookupTable_II[BIASLOOKUPTABLESIZE][4];

public:

  /// Returns the value of the i-th B-spline basis function with segment dependent parameter
  static double B (int, double);

  /// Returns the value of the i-th B-spline basis function with global parameter
  static double N(int i, double u, int L);

  /// Returns the 1st derivative value of the i-th B-spline basis function
  static double B_I(int, double);

  /// Returns the 2nd derivative value of the i-th B-spline basis function
  static double B_II(int, double);

  /// Constructor
  irtkBSplineBiasField();

  /// Constructor
  irtkBSplineBiasField(const irtkGreyImage &image, double dx, double dy, double dz);

  /// Contructor
  irtkBSplineBiasField(const irtkGreyImage &image, int x, int y, int z, bool bounding_box = false, int padding = -1);

  /// Copy Constructor
  irtkBSplineBiasField(const class irtkBSplineBiasField &);

  /// Destructor
  virtual ~irtkBSplineBiasField();

  /// Returns the size of the B-spline lookup table
  virtual int LUTSize() const;

  /** Approximate displacements: This function takes a set of points and a
      set of displacements and find a FFD which approximates these
      displacements. After approximatation the displacements replaced by
      the residual displacement errors at the points */
  virtual double Approximate(double *, double *, double *, double *, int);

  /// Calculate weighted least square fit of B-spline to data
  virtual void WeightedLeastSquares(double *x1, double *y1, double *z1, double *bias, double *weights, int no);

  /** Interpolates displacements: This function takes a set of displacements
      defined at the control points and finds a FFD which interpolates these
      displacements.
      \param dbias The bias at each control point. */
  virtual void Interpolate(double* dbias);

  /// Subdivide FFD
  virtual void Subdivide();

  /// Calculates the bias (for a point in FFD coordinates) with checks
  virtual double FFD1(double, double, double) const;

  /// Calculates the bias  (for a point in FFD coordinates) without checks
  virtual double FFD2(double, double, double) const;

  /// Calculates the bias
  virtual double Bias(double, double, double);

  /// Calculates the bias
  virtual double Bias2(double, double, double);

  /// Reads FFD from file
  virtual void Read (char *);

  /// Writes FFD to file
  virtual void Write(char *);

  /// Returns a string with the name of the instantiated class
  virtual const char *NameOfClass();

  /// Print info
  virtual void Print();

  /// Find index in matrix for calculating least square fit
  virtual int Ind(int, int, int);
};

inline int irtkBSplineBiasField::Ind(int a, int b, int c)
{
  if (a<0) a=0;
  if (b<0) b=0;
  if (c<0) c=0;
  if (a>_x-1) a=_x-1;
  if (b>_y-1) b=_y-1;
  if (c>_z-1) c=_z-1;

  return a+b*_x+c*_y*_x;
}

inline double irtkBSplineBiasField::B(int i, double t)
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

inline double irtkBSplineBiasField::N(int i, double u, int L)
{
  if ((i<0)||(i>L-1)) return 0;
  int l= (int) floor(u);
  if (l==L-1) l=L-2;

  double value=B(i-l+1,u-l);
  if (l==0) {
    if (i==0) value += 2*B(0,u-l);
    if (i==1) value -= B(0,u-l);
  }

  if (l==L-2) {
    if (i==L-1) value += 2*B(3,u-l);
    if (i==L-2) value -= B(3,u-l);
  }

  return value;

}


inline double irtkBSplineBiasField::B_I(int i, double t)
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

inline double irtkBSplineBiasField::B_II(int i, double t)
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

inline double irtkBSplineBiasField::B0(double t)
{
  return (1-t)*(1-t)*(1-t)/6.0;
}

inline double irtkBSplineBiasField::B1(double t)
{
  return (3*t*t*t - 6*t*t + 4)/6.0;
}

inline double irtkBSplineBiasField::B2(double t)
{
  return (-3*t*t*t + 3*t*t + 3*t + 1)/6.0;
}

inline double irtkBSplineBiasField::B3(double t)
{
  return (t*t*t)/6.0;
}

inline double irtkBSplineBiasField::B0_I(double t)
{
  return -(1-t)*(1-t)/2.0;
}

inline double irtkBSplineBiasField::B1_I(double t)
{
  return (9*t*t - 12*t)/6.0;
}

inline double irtkBSplineBiasField::B2_I(double t)
{
  return (-9*t*t + 6*t + 3)/6.0;
}

inline double irtkBSplineBiasField::B3_I(double t)
{
  return (t*t)/2.0;
}

inline double irtkBSplineBiasField::B0_II(double t)
{
  return 1 - t;
}

inline double irtkBSplineBiasField::B1_II(double t)
{
  return 3*t - 2;
}

inline double irtkBSplineBiasField::B2_II(double t)
{
  return -3*t + 1;
}

inline double irtkBSplineBiasField::B3_II(double t)
{
  return t;
}

inline int irtkBSplineBiasField::LUTSize() const
{
  return BIASLOOKUPTABLESIZE;
}

inline double irtkBSplineBiasField::Bias(double x, double y, double z)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  return this->FFD1(u, v, w);
}

inline double irtkBSplineBiasField::Bias2(double x, double y, double z)
{
  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  return this->FFD2(u, v, w);
}

inline const char *irtkBSplineBiasField::NameOfClass()
{
  return "irtkBSplineBiasField";
}

#endif
