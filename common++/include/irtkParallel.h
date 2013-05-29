/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKPARALLEL_H

#define _IRTKPARALLEL_H

/// Debugging level of parallel code
extern int tbb_debug;

/// Number of threads to use
extern int tbb_no_threads;

// If TBB is available and BUILD_TBB_EXE is set to ON, use TBB to execute
// any parallelizable code concurrently
//
// Attention: DO NOT define TBB_DEPRECATED by default or before including the
//            other TBB header files, in particular parallel_for. The deprecated
//            behavior of parallel_for is to not choose the chunk size (grainsize)
//            automatically!
//
// http://software.intel.com/sites/products/documentation/doclib/tbb_sa/help/tbb_userguide/Automatic_Chunking.htm
#ifdef HAS_TBB
#  include <tbb/task_scheduler_init.h>
#  include <tbb/blocked_range.h>
#  include <tbb/blocked_range2d.h>
#  include <tbb/blocked_range3d.h>
#  include <tbb/parallel_for.h>
#  include <tbb/parallel_reduce.h>
#  include <tbb/tick_count.h>
#  include <tbb/concurrent_queue.h>
#  include <tbb/mutex.h>
using namespace tbb;
// Otherwise, use dummy implementations of TBB classes/functions which allows
// developers to write parallelizable code as if TBB was available and yet
// executes the code serially due to the lack of TBB (or BUILD_TBB_EXE set to OFF).
// This avoids code duplication and unnecessary conditional code compilation.
#else

  class task_scheduler_init
  {
  public:
    task_scheduler_init(int) {}
    void terminate() {}
  };

  template <typename T>
  class blocked_range
  {
    int _lbound;
    int _ubound;
  public:
    blocked_range(int l, int u, int = 1) : _lbound(l), _ubound(u) {}
    int begin() const { return _lbound; }
    int end()   const { return _ubound; }
  };

  template <typename T>
  class blocked_range2d
  {
    blocked_range<int> _rows;
    blocked_range<int> _cols;

  public:

    blocked_range2d(int rl, int ru,
                    int cl, int cu)
    :
      _rows (rl, ru),
      _cols (cl, cu)
    {
    }

    blocked_range2d(int rl, int ru, int,
                    int cl, int cu, int)
    :
      _rows (rl, ru),
      _cols (cl, cu)
    {
    }

    const blocked_range<int> &rows() const { return _rows; }
    const blocked_range<int> &cols() const { return _cols; }
  };

  template <typename T>
  class blocked_range3d
  {
    blocked_range<int> _pages;
    blocked_range<int> _rows;
    blocked_range<int> _cols;

  public:

    blocked_range3d(int pl, int pu,
                    int rl, int ru,
                    int cl, int cu)
    :
      _pages(pl, pu),
      _rows (rl, ru),
      _cols (cl, cu)
    {
    }

    blocked_range3d(int pl, int pu, int,
                    int rl, int ru, int,
                    int cl, int cu, int)
    :
      _pages(pl, pu),
      _rows (rl, ru),
      _cols (cl, cu)
    {
    }

    const blocked_range<int> &pages() const { return _pages; }
    const blocked_range<int> &rows() const { return _rows; }
    const blocked_range<int> &cols() const { return _cols; }
  };

  template <class Range, class Body>
  void parallel_for(const Range &range, const Body &body) {
    body(range);
  }

  template <class Range, class Body>
  void parallel_reduce(const Range &range, const Body &body) {
    body(range);
  }
#endif

/// Preprocessor flag to over all remove timing code from binary must be
/// defined non-zero before any include statement if timing should be used
#ifndef USE_TIMING
#  define USE_TIMING 0
#endif

/// Global (non-thread safe) flag to enable/disable timing at run time.
///
/// Should be set in the main function before any processing starts, e.g.,
/// depending on a command-line flag (for example -v -verbose).
/// If less or equal to zero, no timing measurements are printed to screen.
/// Otherwise, whether a timing measure is output or not depends on the
/// set debugging level.
extern int debug_time;

/// Start measurement of execution time of current code block
///
/// @code
/// {
///   IRTK_START_TIMING();
///   // do some work here
///   IRTK_END_TIMING("example section");
/// }
/// @endcode
///
/// @sa IRTK_END_TIMING
#if USE_TIMING && defined(HAS_TBB)
#  define IRTK_START_TIMING()   tick_count t_start = tick_count::now()
#elif USE_TIMING
#  define IRTK_START_TIMING()   clock_t t_start = clock()
#else
#  define IRTK_START_TIMING()   do {} while (false)
#endif

/// End measurement of execution time of current code block
///
/// @code
/// {
///   IRTK_START_TIMING();
///   // do some work here
///   IRTK_END_TIMING("example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at compile time depending on the
///       USE_TIMING flag.
///
/// @sa IRTK_START_TIMING
#if USE_TIMING && defined(HAS_TBB)
#  define IRTK_END_TIMING(section_name) \
     ::std::cout << "CPU time for " section_name ": " \
                 << (tick_count::now() - t_start).seconds() \
                 << " secs" << endl
#elif USE_TIMING
#  define IRTK_END_TIMING(section_name) \
   ::std::cout << "CPU time for " section_name ": " \
               << (static_cast<double>(clock() - t_start) / CLOCKS_PER_SEC) \
               << " secs" << endl
#else
#  define IRTK_END_TIMING(section)   do {} while (false)
#endif

/// End measurement of execution time of current code block
///
/// In the following example, the timing measurement are only printed if
/// the global execution time debugging level is greater or equal the
/// specified debugging level. Hence, if debug_time is 1, only the time of
/// the following example section is printed. If debug_time is greater than 1,
/// also the times needed for each execution of inner block are printed together
/// with the summarizing total execution time of the example section. Otherwise,
/// if debug_time is less than 1, no time measurements are printed at all.
/// @code
/// {
///   IRTK_START_TIMING();
///   // do some work here
///     {
///       IRTK_START_TIMING();
///       // do some part of the work here
///       IRTK_DEBUG_TIMING(2, "example part");
///     }
///   // possibly combine results and maybe do some clean up work here
///   IRTK_DEBUG_TIMING(1, "example section");
/// }
/// @endcode
///
/// @note Whether or not the execution time is actually being measured and
///       printed to screen is decided at runtime depending on the global
///       variable debug_time.
///
/// @sa IRTK_START_TIMING
#if USE_TIMING && defined(HAS_TBB)
#  define IRTK_DEBUG_TIMING(level, section_name) \
     if (debug_time >= level) ::std::cout << "CPU time for " section_name ": " \
                                          << (tick_count::now() - t_start).seconds() \
                                          << " secs" << endl
#elif USE_TIMING
#  define IRTK_DEBUG_TIMING(level, section_name) \
     if (debug_time >= level) ::std::cout << "CPU time for " section_name ": " \
               << (static_cast<double>(clock() - t_start) / CLOCKS_PER_SEC) \
               << " secs" << endl
#else
#  define IRTK_DEBUG_TIMING(level, section)   do {} while (false)
#endif


#endif
