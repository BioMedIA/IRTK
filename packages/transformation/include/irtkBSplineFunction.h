/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKBSPLINEFUNCTION_H

#define _IRTKBSPLINEFUNCTION_H


#include <irtkCommon.h> // round()


/**
 * Cubic B-spline, its basis functions, and their derivatives.
 *
 * This static class provides methods to evaluate the cubic B-spline function,
 * each of its four basis functions, and the corresponding first and second
 * derivatives in the interval [0, 1]. It is in particular intended for use
 * by classes implementing the irtkFreeFormTransformation interface using
 * a cubic B-spline basis to represent the displacement or velocity field,
 * respectively.
 *
 * \note Though not explicitly enforced by protecting these members,
 *       do <b>not</b> modify the lookup tables of precomputed function
 *       values as it will affect all code which makes use of them!
 *       Moreover, make sure to call the static Initialize() method at least
 *       once in your program before accessing any of these precomputed values.
 *
 * Usage:
 * \code
 * irtkBSplineFunction kernel;
 * // evaluate function explicitly
 * double exact = kernel.B1(0.2);
 * // or use lookup table; note that the Initialize() method must only be
 * // called once during the life-time of the program
 * irtkFreeFormTransformationBSplineFunction::Initialize();
 * double approx = kernel.LookupTable[kernel.VariableToIndex(0.2)][1];
 * \endcode
 */

class irtkBSplineFunction
{
public:

  /// Returns the value of the B-spline function
  static double B(double);

  /// Returns the value of the i-th B-spline basis function
  static double B(int, double);

  /// Returns the value of the first B-spline basis function
  static double B0(double);

  /// Returns the value of the second B-spline basis function
  static double B1(double);

  /// Returns the value of the third B-spline basis function
  static double B2(double);

  /// Returns the value of the fourth B-spline basis function
  static double B3(double);

  /// Returns the 1st derivative value of the i-th B-spline basis function
  static double B_I(int, double);
  
  /// Returns the 1st derivative value of the first B-spline basis function
  static double B0_I(double);

  /// Returns the 1st derivative value of the second B-spline basis function
  static double B1_I(double);

  /// Returns the 1st derivative value of the third B-spline basis function
  static double B2_I(double);

  /// Returns the 1st derivative value of the fourth B-spline basis function
  static double B3_I(double);

  /// Returns the 2nd derivative value of the i-th B-spline basis function
  static double B_II(int, double);

  /// Returns the 2nd derivative value of the first B-spline basis function
  static double B0_II(double);

  /// Returns the 2nd derivative value of the second B-spline basis function
  static double B1_II(double);

  /// Returns the 2nd derivative value of the third B-spline basis function
  static double B2_II(double);

  /// Returns the 2nd derivative value of the fourth B-spline basis function
  static double B3_II(double);

  /// Size of lookup tables used to store pre-computed function values
  static const int LookupTableSize = 1000;

  /// Initialize lookup tables
  static void Initialize(bool = false);

  /// Returns the lookup table index for a given value in [0,1]
  static int VariableToIndex(double);

  /// Lookup table of B-spline basis function values
  static double LookupTable[LookupTableSize][4];

  /// Lookup table of B-spline basis function 1st derivative values
  static double LookupTable_I[LookupTableSize][4];

  /// Lookup table of B-spline basis function 2nd derivative values
  static double LookupTable_II[LookupTableSize][4];

protected:
  
  /// Flag which indicates whether the lookup tables are initialized
  static bool _initialized;
};

inline double irtkBSplineFunction::B(double t)
{
  double value = 0.0;
  t = fabs(t);
  if (t < 2.0) {
    if (t < 1.0) {
      value = 2.0/3.0 + (0.5*t-1.0)*t*t;
    } else {
      t -= 2.0;
      value = -t*t*t/6.0;
    }
  }
  return value;
}


inline double irtkBSplineFunction::B(int i, double t)
{
  switch (i) {
    case 0:  return B0(t);
    case 1:  return B1(t);
    case 2:  return B2(t);
    case 3:  return B3(t);
    default: return 0.0;
  }
}

inline double irtkBSplineFunction::B0(double t)
{
  return (1.0-t)*(1.0-t)*(1.0-t)/6.0;
}

inline double irtkBSplineFunction::B1(double t)
{
  return (3.0*t*t*t - 6.0*t*t + 4.0)/6.0;
}

inline double irtkBSplineFunction::B2(double t)
{
  return (-3.0*t*t*t + 3.0*t*t + 3.0*t + 1.0)/6.0;
}

inline double irtkBSplineFunction::B3(double t)
{
  return (t*t*t)/6.0;
}


inline double irtkBSplineFunction::B_I(int i, double t)
{
  switch (i) {
    case 0:  return B0_I(t);
    case 1:  return B1_I(t);
    case 2:  return B2_I(t);
    case 3:  return B3_I(t);
    default: return 0.0;
  }
}

inline double irtkBSplineFunction::B0_I(double t)
{
  return -(1.0-t)*(1.0-t)/2.0;
}

inline double irtkBSplineFunction::B1_I(double t)
{
  return (9.0*t*t - 12.0*t)/6.0;
}

inline double irtkBSplineFunction::B2_I(double t)
{
  return (-9.0*t*t + 6.0*t + 3.0)/6.0;
}

inline double irtkBSplineFunction::B3_I(double t)
{
  return (t*t)/2.0;
}


inline double irtkBSplineFunction::B_II(int i, double t)
{
  switch (i) {
    case 0:  return B0_II(t);
    case 1:  return B1_II(t);
    case 2:  return B2_II(t);
    case 3:  return B3_II(t);
    default: return 0.0;
  }
}

inline double irtkBSplineFunction::B0_II(double t)
{
  return 1.0 - t;
}

inline double irtkBSplineFunction::B1_II(double t)
{
  return 3.0*t - 2.0;
}

inline double irtkBSplineFunction::B2_II(double t)
{
  return -3.0*t + 1.0;
}

inline double irtkBSplineFunction::B3_II(double t)
{
  return t;
}


inline void irtkBSplineFunction::Initialize(bool force)
{
  if (!_initialized || force) {
    for (int i = 0; i < LookupTableSize; i++) {
      for (int l = 0; l < 4; l++) {
        LookupTable   [i][l] = B   (l, static_cast<double>(i) / static_cast<double>(LookupTableSize-1));
        LookupTable_I [i][l] = B_I (l, static_cast<double>(i) / static_cast<double>(LookupTableSize-1));
        LookupTable_II[i][l] = B_II(l, static_cast<double>(i) / static_cast<double>(LookupTableSize-1));
      }
    }
    _initialized = true;
  }
}

inline int irtkBSplineFunction::VariableToIndex(double t)
{
  return round(t * static_cast<double>(LookupTableSize-1));
}


#endif
