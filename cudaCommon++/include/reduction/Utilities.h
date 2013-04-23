/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : 
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : 
  Version   : 
  Changes   : $Author: bkainz $
  GPL licensed file: Original reduction implementation by: Noel Lopes, GPUMLib

=========================================================================*/

#ifndef GPUMLib_Utilities_h
#define GPUMLib_Utilities_h

#include "CudaDefinitions.h"
#include "../memory/HostArray.h"

//! \addtogroup commonframework Common framework
//! @{

//! Finds the number of threads (multiple of 2) per block that either is greater that the number of threads needed or identical to the maximum number of threads per block.
//! \param threads Number of threads.
//! \param maxThreadsPerBlock Maximum number of threads.
//! \return The number of threads (multiple of 2) per block that either is greater that the number of threads needed or identical to the maximum number of threads per block.
//! \sa MAX_THREADS_PER_BLOCK, NumberBlocks
inline int NumberThreadsPerBlockThatBestFit(int threads, int maxThreadsPerBlock = MAX_THREADS_PER_BLOCK) {
	int nt = 1;
	while(nt < threads && nt < maxThreadsPerBlock) nt <<= 1;

	return nt;
}

//! Finds the number of blocks needed to execute the number of threads specified, given a block size.
//! \param threads Number of threads.
//! \param blockSize Block size.
//! \return The number of blocks needed to execute the number of threads specified.
//! \sa NumberThreadsPerBlockThatBestFit, MAX_THREADS_PER_BLOCK
inline int NumberBlocks(int threads, int blockSize) {
	int nb = threads / blockSize;

	if (threads % blockSize != 0) nb++;

	return nb;
}

//! @}

#endif