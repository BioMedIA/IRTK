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

#ifndef _IRTKGAUSSIAN_H

#define _IRTKGAUSSIAN_H

#include <irtkImage.h>

/**

multivariate gaussian probability distribution

*/

class irtkGaussian : public irtkObject
{

protected:

  double _mi;
  double _sigma;
  double _norm;

public:

  void Initialise(const double &mi, const double &sigma);

  double Evaluate(const double &x);

  GetMacro(norm,double);
};

inline double irtkGaussian::Evaluate(const double &x)
{
  return _norm * exp(-((x - _mi) * (x - _mi)) / (2.0 * _sigma));
}

#endif
