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

void swap16(char *a, char *b, long n)
{
  int i;
  char c;

  for (i = 0; i < n * 2; i += 2) {
    c = a[i];
    a[i] = b[i+1];
    b[i+1] = c;
  }
}

void swap32(char *a, char *b, long n)
{
  int i;
  char c;

  for (i = 0; i < n * 4; i += 4) {
    c = a[i];
    a[i] = b[i+3];
    b[i+3] = c;
    c = a[i+1];
    a[i+1] = b[i+2];
    b[i+2] = c;
  }
}

void swap64(char *a, char *b, long n)
{
  int i;
  char c;

  for (i = 0; i < n * 8; i += 8) {
    c = a[i];
    a[i] = b[i+7];
    b[i+7] = c;
    c = a[i+1];
    a[i+1] = b[i+6];
    b[i+6] = c;
    c = a[i+2];
    a[i+2] = b[i+5];
    b[i+5] = c;
    c = a[i+3];
    a[i+3] = b[i+4];
    b[i+4] = c;
  }
}
