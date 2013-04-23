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

#ifndef GPUMLib_CudaStreams_h
#define GPUMLib_CudaStreams_h

#include <cuda_runtime.h>

namespace GPUMLib {

//! \addtogroup commonframework Common framework
//! @{

//! Represents a CUDA stream
class CudaStream {
	private:
		cudaStream_t stream;

	public:
		//! Constructor
		CudaStream() {
			cudaStreamCreate(&stream);
		}

		//! Destructor
		~CudaStream() {
			cudaStreamDestroy(stream);
		}

		//! \return Returns the stream as a \b cudaStream_t type
		operator cudaStream_t () {
			return stream;
		}
};

//! @}

}

#endif