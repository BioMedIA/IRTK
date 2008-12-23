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
