/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKBSPLINEINTERPOLATEIMAGEFUNCTION_H

#define _IRTKBSPLINEINTERPOLATEIMAGEFUNCTION_H

/**
 * Class for B-spline interpolation of images
 *
 * This class defines and implements the B-spline interpolation of images.
 * Currently supports B-splines of degree 2 to degree 5. By default splines
 * of degree 3 (cubic) are used. For more details see:
 *
 * M. Unser, "Splines: A Perfect Fit for Signal and Image Processing," IEEE
 * Signal Processing Magazine, vol. 16, no. 6, pp. 22-38, November 1999.
 *
 */

class irtkBSplineInterpolateImageFunction : public irtkInterpolateImageFunction
{

private:

  /// Degree of spline
  int _SplineDegree;

  /// Dimension of input image in X-direction
  int _x;

  /// Dimension of input image in Y-direction
  int _y;

  /// Dimension of input image in Z-direction
  int _z;

  /// Dimension of input image in T-direction
  int _t;

  /// Dimension of input image in X-direction minus spline degree
  int _xhalf;

  /// Dimension of input image in Y-direction minus spline degree
  int _yhalf;

  /// Dimension of input image in Z-direction minus spline degree
  int _zhalf;

  /// Image of spline coefficient
  irtkRealImage _coeff;

  /// Initialize anti-causal coefficients
  static double InitialAntiCausalCoefficient(double *, int, double z);

  /// Initialize causal coefficients
  static double InitialCausalCoefficient(double *, int, double, double);

  /// Convert voxel values to B-spline coefficients
  static void ConvertToInterpolationCoefficients(double *, int, double *, int, double);

  /// Compute B-spline coefficients
  void ComputeCoefficients();

public:

  /// Constructor (default spline degree is cubic B-splines)
  irtkBSplineInterpolateImageFunction(int = 3);

  /// Destructor
  ~irtkBSplineInterpolateImageFunction();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Initialize
  virtual void Initialize();

  /// Sets the spline degree (supported are spline degrees 2 to 5)
  virtual void PutSplineDegree(int SplineDegree);

  /// Gets the spline degree
  virtual int GetSplineDegree();

  /// Evaluate the filter at an arbitrary image location (in pixels)
  virtual double Evaluate(double, double, double, double = 0);

  /** Evaluate the filter at an arbitrary image location (in pixels) without
   *  handling boundary conditions. This version is faster than the method
   *  above, but is only defined inside the image domain. */
  virtual double EvaluateInside(double, double, double, double = 0);

};

inline void irtkBSplineInterpolateImageFunction::PutSplineDegree(int SplineDegree)
{
  if ((SplineDegree < 2) || (SplineDegree > 5)) {
    cerr << "irtkBSplineInterpolateImageFunction::PutSplineDegree: Unsupported spline degree\n";
    exit(1);
  }
  _SplineDegree = SplineDegree;
}

inline int irtkBSplineInterpolateImageFunction::GetSplineDegree()
{
  return _SplineDegree;
}

#endif
