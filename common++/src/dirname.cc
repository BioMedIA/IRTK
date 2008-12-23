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

const char *dirname2(const char *path)
{
  static const char dot[] = ".";
  char *last_slash;

  /* Find last '/'.  */
  last_slash = path != NULL ? strrchr (const_cast<char*>(path), '/') : NULL;

  if (last_slash == path)
    /* The last slash is the first character in the string.  We have to
       return "/".  */
    ++last_slash;
  else if (last_slash != NULL && last_slash[1] == '\0')
    /* The '/' is the last character, we have to look further.  */
    last_slash = (char *)memchr(path, last_slash - path, '/');

  if (last_slash != NULL)
    /* Terminate the path.  */
    last_slash[0] = '\0';
  else
    /* This assignment is ill-designed but the XPG specs require to
       return a string containing "." in any case no directory part is
       found and so a static and constant string is required.  */
    path = (char *) dot;

  return path;
}
