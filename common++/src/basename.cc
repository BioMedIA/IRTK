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

#include <irtkCommon.h>

// Copied from glibc-2.1.1

char *basename2(const char *filename)
{
  char *p = strrchr (const_cast<char*>(filename), '/');
  return p ? p + 1 : (char *) filename;
}
