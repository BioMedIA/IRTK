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

#include <irtkGaussian.h>

#include <irtkGeometry.h>

void irtkGaussian::Initialise(const double &mi, const double &sigma)
{
  _mi = mi;
  _sigma = sigma;
  _norm = 1.0 / (sqrt(sigma) * sqrt(M_PI*2.0));
}


