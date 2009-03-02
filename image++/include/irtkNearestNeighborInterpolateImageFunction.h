/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKNEARESTNEIGHBORINTERPOLATEIMAGEFUNCTION_H

#define _IRTKNEARESTNEIGHBORINTERPOLATEIMAGEFUNCTION_H

/**
 * Class for nearest neighbor interpolation of images
 *
 * This class defines and implements the nearest neighbor interpolation  of
 * images.
 */

class irtkNearestNeighborInterpolateImageFunction : public irtkInterpolateImageFunction
{

private:

  /// Dimension of input image in X-direction
  int _x;

  /// Dimension of input image in Y-direction
  int _y;

  /// Dimension of input image in Z-direction
  int _z;

public:

  /// Constructor
  irtkNearestNeighborInterpolateImageFunction();

  /// Destructor
  ~irtkNearestNeighborInterpolateImageFunction();

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
