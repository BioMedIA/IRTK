/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkGaussian.h>

#include <irtkGeometry.h>

void irtkGaussian::Initialise(const double &mi, const double &sigma)
{
  _mi = mi;
  _sigma = sigma;
  _norm = 1.0 / (sqrt(sigma) * sqrt(M_PI*2.0));
}


