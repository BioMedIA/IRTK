/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSINCINTERPOLATEIMAGEFUNCTION2D_H

#define _IRTKSINCINTERPOLATEIMAGEFUNCTION2D_H

/**
 * Class for sinc interpolation of images
 *
 * This class defines and implements the sinc interpolation of
 * images.
 */

class irtkSincInterpolateImageFunction2D : public irtkInterpolateImageFunction
{

private:

  /// Dimension of input image in X-direction
  int _x;

  /// Dimension of input image in Y-direction
  int _y;

public:

  /// Constructor
  irtkSincInterpolateImageFunction2D();

  /// Destructor
  ~irtkSincInterpolateImageFunction2D();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Initialize
  virtual void Initialize();

  /// Evaluate the filter at an arbitrary image location (in pixels)
  virtual double Evaluate(double, double, double = 0, double = 0);

  /** Evaluate the filter at an arbitrary image location (in pixels) without
   *  handling boundary conditions. This version is faster than the method
   *  above, but is only defined inside the image domain. */
  virtual double EvaluateInside(double, double, double = 0, double = 0);

};

#endif
