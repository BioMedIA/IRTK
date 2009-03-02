/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKLINEARINTERPOLATEIMAGEFUNCTION_H

#define _IRTKLINEARINTERPOLATEIMAGEFUNCTION_H

/**
 * Class for linear interpolation of images
 *
 * This class defines and implements the linear interpolation of images.
 *
 */

class irtkLinearInterpolateImageFunction : public irtkInterpolateImageFunction
{

private:

  /// Dimension of input image in X-direction
  int _x;

  /// Dimension of input image in Y-direction
  int _y;

  /// Dimension of input image in Z-direction
  int _z;

  /// Offsets for fast pixel access
  int _offset1, _offset2, _offset3, _offset4;
  int _offset5, _offset6, _offset7, _offset8;

public:

  /// Constructor
  irtkLinearInterpolateImageFunction();

  /// Destructor
  ~irtkLinearInterpolateImageFunction();

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
