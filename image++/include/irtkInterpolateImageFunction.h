/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKINTERPOLATEIMAGEFUNCTION_H

#define _IRTKINTERPOLATEIMAGEFUNCTION_H

/// Image interpolation modes
typedef enum { Interpolation_NN,
               Interpolation_Linear,
               Interpolation_BSpline,
               Interpolation_CSpline,
               Interpolation_Sinc,
               Interpolation_Gaussian
             } irtkInterpolationMode;

/**
 * Abstract base class for any general image interpolation filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and sample that image at arbitrary
 * location. Each derived class has to implement all abstract member functions.
 */

template <class VoxelType> class irtkInterpolateImageFunction : public irtkImageFunction<VoxelType>
{

protected:

  /// Flag whether to clamp interpolated values to min and max range of
  /// original data
  int _clamped;

  /// Min and max range of original data
  VoxelType _min, _max;

  /** Domain of input image for which this image interpolation function can
   *  be used without handling any form of boundary conditions. */
  double _x1, _y1, _z1, _x2, _y2, _z2;

public:

  /// Constructor
  irtkInterpolateImageFunction();

  /// Deconstuctor
  virtual ~irtkInterpolateImageFunction();

  /// Constuctor
  static irtkInterpolateImageFunction *New(irtkInterpolationMode,
      irtkBaseImage *);

  /// Returns the name of the class
  /// Initialize
  virtual void Initialize();

  /// Evaluate the filter at an arbitrary image location (in pixels)
  virtual double Evaluate(double, double, double, double = 0) = 0;

  /** Evaluate the filter at an arbitrary image location (in pixels) without
   *  handling boundary conditions. This version is faster than the method
   *  above, but is only defined inside the image domain. */
  virtual double EvaluateInside(double, double, double, double = 0) = 0;

  /** Check if the location is inside the image domain for which this image
   *  interpolation function can be used without handling any form of boundary
   *  conditions. */
  virtual Bool IsInside(double, double, double);

  /** Returns the image domain for which this image interpolation function
      can be used without handling any form of boundary conditions. */
  virtual void Inside(double &, double &, double &,
                      double &, double &, double &);

};

template <class VoxelType> inline Bool irtkInterpolateImageFunction<VoxelType>::IsInside(double x, double y, double z)
{
  return ((x > _x1) && (x < _x2) && (y > _y1) && (y < _y2) &&
          (z > _z1) && (z < _z2));
}

template <class VoxelType> inline void irtkInterpolateImageFunction<VoxelType>::Inside(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2)
{
  x1 = _x1;
  y1 = _y1;
  z1 = _z1;
  x2 = _x2;
  y2 = _y2;
  z2 = _z2;
}

#include <irtkNearestNeighborInterpolateImageFunction.h>
#include <irtkLinearInterpolateImageFunction.h>
#include <irtkBSplineInterpolateImageFunction.h>
#include <irtkCSplineInterpolateImageFunction.h>
#include <irtkSincInterpolateImageFunction.h>
#include <irtkGaussianInterpolateImageFunction.h>

#include <irtkNearestNeighborInterpolateImageFunction2D.h>
#include <irtkLinearInterpolateImageFunction2D.h>
#include <irtkBSplineInterpolateImageFunction2D.h>
#include <irtkCSplineInterpolateImageFunction2D.h>
#include <irtkSincInterpolateImageFunction2D.h>
#include <irtkGaussianInterpolateImageFunction2D.h>

#endif
