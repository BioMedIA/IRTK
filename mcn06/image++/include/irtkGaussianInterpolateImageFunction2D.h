/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGAUSSIANINTERPOLATEIMAGEFUNCTION2D_H

#define _IRTKGAUSSIANINTERPOLATEIMAGEFUNCTION2D_H

/**
 * Class for linear interpolation of images
 *
 * This class defines and implements the linear interpolation of images.
 *
 */

class irtkGaussianInterpolateImageFunction2D : public irtkInterpolateImageFunction
{

private:

  /// Dimension of input image in X-direction
  int _x;

  /// Dimension of input image in Y-direction
  int _y;

  /// Dimension of input image in Z-direction
  int _z;

  /// Dimension of voxels in input image in X-direction
  double _xsize;

  /// Dimension of voxels in input image in Y-direction
  double _ysize;

  /// Dimension of voxels in input image in Z-direction
  double _zsize;

  /// Filter size in X-direction
  double _ExtentX;

  /// Filter size in Y-direction
  double _ExtentY;

  /// Sigma (in mm)
  double _Sigma;

public:

  /// Constructor
  irtkGaussianInterpolateImageFunction2D(double sigma = 1);

  /// Destructor
  ~irtkGaussianInterpolateImageFunction2D();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Initialize
  virtual void Initialize();

  /// Evaluate
  virtual double Evaluate(double, double, double, double = 0);

  /** Evaluate the filter at an arbitrary image location (in pixels) without
   *  handling boundary conditions. This version is faster than the method
   *  above, but is only defined inside the image domain. */
  virtual double EvaluateInside(double, double, double, double = 0);

};

#endif

