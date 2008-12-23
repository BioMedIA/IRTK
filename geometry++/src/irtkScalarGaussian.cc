/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkGeometry.h>

irtkScalarGaussian::irtkScalarGaussian()
{
  _Sigma_x = 1;
  _Sigma_y = 1;
  _Sigma_z = 1;
  _X_0 = 0;
  _Y_0 = 0;
  _Z_0 = 0;
  _Norm = 1.0 / (sqrt(2.0 * M_PI) * _Sigma_x *
                 sqrt(2.0 * M_PI) * _Sigma_y *
                 sqrt(2.0 * M_PI) * _Sigma_z);
}

irtkScalarGaussian::irtkScalarGaussian(double _sigma)
{
  _Sigma_x = _sigma;
  _Sigma_y = _sigma;
  _Sigma_z = _sigma;
  _X_0 = 0;
  _Y_0 = 0;
  _Z_0 = 0;
  _Norm = 1.0 / (sqrt(2.0 * M_PI) * _Sigma_x *
                 sqrt(2.0 * M_PI) * _Sigma_y *
                 sqrt(2.0 * M_PI) * _Sigma_z);
}

irtkScalarGaussian::irtkScalarGaussian(double _sigma, double x_0, double y_0, double z_0)
{
  _Sigma_x = _sigma;
  _Sigma_y = _sigma;
  _Sigma_z = _sigma;
  _X_0 = x_0;
  _Y_0 = y_0;
  _Z_0 = z_0;
  _Norm = 1.0 / (sqrt(2.0 * M_PI) * _Sigma_x *
                 sqrt(2.0 * M_PI) * _Sigma_y *
                 sqrt(2.0 * M_PI) * _Sigma_z);
}

irtkScalarGaussian::irtkScalarGaussian(double _sigma_x, double _sigma_y, double _sigma_z, double x_0, double y_0, double z_0)
{
  _Sigma_x = _sigma_x;
  _Sigma_y = _sigma_y;
  _Sigma_z = _sigma_z;
  _X_0 = x_0;
  _Y_0 = y_0;
  _Z_0 = z_0;
  if ((_Sigma_x != 0) && (_Sigma_y != 0) && (_Sigma_z != 0)) {
    _Norm = 1.0 / (sqrt(2.0 * M_PI) * _Sigma_x *
                   sqrt(2.0 * M_PI) * _Sigma_y *
                   sqrt(2.0 * M_PI) * _Sigma_z);
  } else {
    cerr << "irtkScalarGaussian::irtkScalarGaussian: Warning, divide by zero" << endl;
    _Norm = 0.0;
  }
}

// Destructor
irtkScalarGaussian::~irtkScalarGaussian()
{
}

// Evaluate
double irtkScalarGaussian::Evaluate(double x, double y, double z)
{
  return _Norm * exp(- ((x-_X_0) * (x-_X_0))/(2.0 * _Sigma_x * _Sigma_x)
                     - ((y-_Y_0) * (y-_Y_0))/(2.0 * _Sigma_y * _Sigma_y)
                     - ((z-_Z_0) * (z-_Z_0))/(2.0 * _Sigma_z * _Sigma_z));
}
