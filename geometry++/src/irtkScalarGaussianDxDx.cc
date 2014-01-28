/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkGeometry.h>

#include <irtkScalarGaussian.h>
#include <irtkScalarGaussianDxDx.h>

irtkScalarGaussianDxDx::irtkScalarGaussianDxDx(double _sigma_x, double _sigma_y, double _sigma_z, double x_0, double y_0, double z_0)
{
  _Sigma_x = _sigma_x;
  _Sigma_y = _sigma_y;
  _Sigma_z = _sigma_z;
  _X_0 = x_0;
  _Y_0 = y_0;
  _Z_0 = z_0;
  _VarX = _Sigma_x * _Sigma_x;
  _VarY = _Sigma_y * _Sigma_y;
  _VarZ = _Sigma_z * _Sigma_z;

  if ((_Sigma_x != 0) && (_Sigma_y != 0) && (_Sigma_z != 0)){
    _Factor = 1.0 / (sqrt(2.0 * M_PI) * _Sigma_x * 
		   sqrt(2.0 * M_PI) * _Sigma_y * 
		   sqrt(2.0 * M_PI) * _Sigma_z);
  } else {
    cerr << "irtkScalarGaussianDxDx::irtkScalarGaussianDxDx: Warning, divide by zero" << endl;
    _Factor = 0.0;
  }
}

irtkScalarGaussianDxDx::~irtkScalarGaussianDxDx()
{

}

double irtkScalarGaussianDxDx::Evaluate(double x, double y, double z)
{
  _Exp = exp(- ((x-_X_0) * (x-_X_0))/(2.0 * _VarX) 
		     - ((y-_Y_0) * (y-_Y_0))/(2.0 * _VarY) 
		     - ((z-_Z_0) * (z-_Z_0))/(2.0 * _VarZ));
  //return _VarX * (- (1.0 / _VarX) * _Factor * _Exp + (x - _X_0) * (1.0 / _VarX) * (x - _X_0) * (1.0 / _VarX) * _Factor * _Exp);
  return - (1.0 / _VarX) * _Factor * _Exp + (x - _X_0) * (1.0 / _VarX) * (x - _X_0) * (1.0 / _VarX) * _Factor * _Exp;
}
