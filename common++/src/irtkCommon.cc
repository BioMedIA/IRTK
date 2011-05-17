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

#ifdef HAS_TBB

// Default: Debugging of TBB code
int tbb_debug = false;

// Default: Number of threads is determined automatically
int tbb_no_threads = task_scheduler_init::automatic;

#endif


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

