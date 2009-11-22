/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkCommon.h>

// Copied from glibc-2.1.1

char *basename2(const char *filename)
{
  char *p = strrchr (const_cast<char*>(filename), '/');
  return p ? p + 1 : (char *) filename;
}
