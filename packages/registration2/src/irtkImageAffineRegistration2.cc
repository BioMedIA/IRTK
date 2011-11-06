/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration2.h>

#include <irtkGradientImageFilter.h>

#include <irtkHomogeneousTransformationIterator.h>

#ifdef WIN32
#include <time.h>
#else
#include <sys/resource.h>
#endif

