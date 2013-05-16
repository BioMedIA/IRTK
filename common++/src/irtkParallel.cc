/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkParallel.h>

// Default: No debugging of execution time
int debug_time = 0;

// Default: No debugging of TBB code
int tbb_debug = false;

// Default: Number of threads is determined automatically
#ifdef HAS_TBB
int tbb_no_threads = task_scheduler_init::automatic;
#else
int tbb_no_threads = 1;
#endif
