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

void PrintVersion(ostream &out, const char* revisionString)
{
	out << "(SVN revision: ";
	// Extract the number from the SVN supplied revision string.
	int len = strlen(revisionString);
	const char *ptr = revisionString;
	for (int i = 0; i < len; ++i){
		if (*ptr >= '0' && *ptr <= '9')
				out << *ptr;
		++ptr;
	}
	out << ")\n";

}

