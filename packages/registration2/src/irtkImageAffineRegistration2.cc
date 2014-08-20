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

#include <irtkRegistration2.h>

#include <irtkGradientImageFilter.h>

#include <irtkHomogeneousTransformationIterator.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

