/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSCALARGAUSSIANDY_H

#define _IRTKSCALARGAUSSIANDY_H

/**

  Class for computing Gaussian first derivative with respect to y
  
*/

class irtkScalarGaussianDy : public irtkScalarGaussian {
  
protected:

  /// Normalization of the Gaussian function
  double _Factor;
  double _Exp;
  double _VarX;
  double _VarY;
  double _VarZ;

public:

  /// Constructor with anisotropic sigma and specific center
  irtkScalarGaussianDy(double sigma_x, double sigma_y, double sigma_z, double x_0, double y_0, double z_0);

  /// Virtual destructor
  virtual ~irtkScalarGaussianDy();

  /// Virtual local evaluation function
  double Evaluate(double, double, double);
  
};

#endif




