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

std::ostream& operator<<(std::ostream& os, const irtkException& ex)
{
  os << "File name: ";
  if (ex.GetFileName() != "")
    os << ex.GetFileName();
  else
    os << "UNKNOWN";
  os << std::endl;

  os << "Line: ";
  if (ex.GetLine() != 0)
    os << ex.GetLine();
  else
    os << "UNKNOWN";
  os << std::endl;

  os << "Message: " << ex.GetMessage();

  return os;
}
