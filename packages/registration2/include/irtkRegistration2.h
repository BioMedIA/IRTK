/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#ifndef _IRTKREGISTRATION2_H

#define _IRTKREGISTRATION2_H

#include <irtkRegistration.h>

class irtkRegistration2 : public irtkObject
{

public:

  /// Evaluate similarity metric
  virtual double Evaluate() = 0;

  /// Evaluate gradient of similarity metric
  virtual double EvaluateGradient(double *) = 0;

};

#include <irtkImageRegistration2.h>
#include <irtkMultipleImageRegistration2.h>

#endif
