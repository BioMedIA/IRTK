/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSCALARGAUSSIANDXDY_H

#define _IRTKSCALARGAUSSIANDXDY_H

/**

  Class for computing Gaussian second derivative with respect to x and y
  
*/

class irtkScalarGaussianDxDy : public irtkScalarGaussian {
  
protected:

  /// Normalization of the Gaussian function
  double _Factor;
  double _Exp;
  double _VarX;
  double _VarY;
  double _VarZ;

public:

  /// Constructor with anisotropic sigma and specific center
  irtkScalarGaussianDxDy(double sigma_x, double sigma_y, double sigma_z, double x_0, double y_0, double z_0);

  /// Virtual destructor
  virtual ~irtkScalarGaussianDxDy();

  /// Virtual local evaluation function
  double Evaluate(double, double, double);
  
};

#endif
