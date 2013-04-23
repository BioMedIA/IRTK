/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : 
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2011 onwards
  Date      : 
  Version   : 
  Changes   : $Author: bkainz $

=========================================================================*/

#define _IMPLEMENTS_CUGENERICIMAGE_

#include <irtkImage.h>

#include <irtkFileToImage.h>
#include <irtkImageToFile.h>

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>

//testing
//#undef USE_CUDA_NPP
#define USE_CUDA_NPP

#ifdef USE_CUDA_NPP
#include <npp.h>
#include "nppExceptions.h"
#endif

#ifdef DUSE_GPL_CODE_ON
#include "reduction/reduction.h"
#endif

#include "irtkCUGenericImage.h"
#include "irtkCUImage_kernels.cuh"

static bool cudaInitDone = false;

template <class VoxelType> irtkCUGenericImage<VoxelType>::irtkCUGenericImage() : irtkGenericImage<VoxelType>()
{
	cudaCheckSets();
	this->_attr._x = 0;
	this->_attr._y = 0;
	this->_attr._z = 0;
	this->_attr._t = 0;

	// Initialize data
	this->_matrix  = NULL;
	d_matrix = NULL;
}


template <class VoxelType> irtkCUGenericImage<VoxelType>::irtkCUGenericImage(int x, int y, int z, int t) : irtkGenericImage<VoxelType>()
{
	cudaCheckSets();
	irtkImageAttributes attr;

	attr._x = x;
	attr._y = y;
	attr._z = z;
	attr._t = t;

	// Initialize data
	this->_matrix = NULL;
	d_matrix = NULL;
	// Initialize rest of class
	this->Initialize(attr);

	//map host memory to device 
	checkCudaErrors(cudaHostRegister(this->GetScalarPointer(),this->GetNumberOfVoxels()*sizeof(VoxelType), cudaHostRegisterMapped));
	cudaDeviceSynchronize();
	checkCudaErrors(cudaHostGetDevicePointer(&d_matrix, this->GetPointerToVoxels(), 0));
	cudaDeviceSynchronize();
}

template <class VoxelType> irtkCUGenericImage<VoxelType>::irtkCUGenericImage(char *filename)
{
	cudaCheckSets();
	// Initialize data
	this->_matrix = NULL;
	d_matrix = NULL;

	// Read image
	this->Read(filename);

	//map host memory to device 
	checkCudaErrors(cudaHostRegister(this->GetPointerToVoxels(),this->GetNumberOfVoxels()*sizeof(VoxelType), cudaHostRegisterMapped));
	cudaDeviceSynchronize();
	checkCudaErrors(cudaHostGetDevicePointer(&d_matrix, this->GetPointerToVoxels(), 0));
	cudaDeviceSynchronize();
}


template <class VoxelType> irtkCUGenericImage<VoxelType>::irtkCUGenericImage(const irtkImageAttributes &attr) : irtkGenericImage<VoxelType>()
{
	cudaCheckSets();
	// Initialize data
	this->_matrix  = NULL;
	d_matrix = NULL;

	// Initialize rest of class
	this->Initialize(attr);

	//map host memory to device 
	checkCudaErrors(cudaHostRegister(this->GetPointerToVoxels(),this->GetNumberOfVoxels()*sizeof(VoxelType), cudaHostRegisterMapped));
	cudaDeviceSynchronize();
	checkCudaErrors(cudaHostGetDevicePointer(&d_matrix, this->GetPointerToVoxels(), 0));
	cudaDeviceSynchronize();
}

template <class VoxelType> irtkCUGenericImage<VoxelType>::irtkCUGenericImage(const irtkCUGenericImage &image) : irtkGenericImage<VoxelType>()
{
	cudaCheckSets();
	int i, n;
	VoxelType *ptr1, *ptr2;

	// Initialize data
	this->_matrix = NULL;
	d_matrix = NULL;
	// Initialize rest of class
	this->Initialize(image._attr);

	n    = this->GetNumberOfVoxels();
	ptr1 = this->GetPointerToVoxels();
	ptr2 = image.GetPointerToVoxels();
	//TODO memcpy
	for (i = 0; i < n; i++) {
		ptr1[i] = ptr2[i];
	}

	//map host memory to device 
	checkCudaErrors(cudaHostRegister(this->GetPointerToVoxels(),this->GetNumberOfVoxels()*sizeof(VoxelType), cudaHostRegisterMapped));
	cudaDeviceSynchronize();
	checkCudaErrors(cudaHostGetDevicePointer(&d_matrix, this->GetPointerToVoxels(), 0));
	cudaDeviceSynchronize();
}

template <class VoxelType> template <class VoxelType2> irtkCUGenericImage<VoxelType>::irtkCUGenericImage(const irtkGenericImage<VoxelType2> &image)
{
	cudaCheckSets();
	int i, n;
	VoxelType  *ptr1;
	VoxelType2 *ptr2;

	// Initialize data
	this->_matrix = NULL;
	d_matrix = NULL;
	// Initialize rest of class
	this->Initialize(image.GetImageAttributes());

	n    = this->GetNumberOfVoxels();
	ptr1 = this->GetPointerToVoxels();
	ptr2 = image.GetPointerToVoxels();
	//TODO memcpy
	for (i = 0; i < n; i++) {
		ptr1[i] = static_cast<VoxelType>(ptr2[i]);
	}

	//map host memory to device 
	checkCudaErrors(cudaHostRegister(this->GetPointerToVoxels(),this->GetNumberOfVoxels()*sizeof(VoxelType), cudaHostRegisterMapped));
	cudaDeviceSynchronize();
	checkCudaErrors(cudaHostGetDevicePointer(&d_matrix, this->GetPointerToVoxels(), 0));
	cudaDeviceSynchronize();
}

template <class VoxelType> template <class VoxelType2> irtkCUGenericImage<VoxelType>::irtkCUGenericImage(const irtkCUGenericImage<VoxelType2> &image)
{
	cudaCheckSets();
	int i, n;
	VoxelType  *ptr1;
	VoxelType2 *ptr2;

	// Initialize data
	this->_matrix = NULL;
	d_matrix = NULL;
	// Initialize rest of class
	this->Initialize(image.GetImageAttributes());

	n    = this->GetNumberOfVoxels();
	ptr1 = this->GetPointerToVoxels();
	ptr2 = image.GetPointerToVoxels();
	//TODO memcpy
	for (i = 0; i < n; i++) {
		ptr1[i] = static_cast<VoxelType>(ptr2[i]);
	}

	//map host memory to device 
	checkCudaErrors(cudaHostRegister(this->GetPointerToVoxels(),this->GetNumberOfVoxels()*sizeof(VoxelType), cudaHostRegisterMapped));
	cudaDeviceSynchronize();
	checkCudaErrors(cudaHostGetDevicePointer(&d_matrix, this->GetPointerToVoxels(), 0));
	cudaDeviceSynchronize();
}

template <class VoxelType> irtkCUGenericImage<VoxelType>::~irtkCUGenericImage(void)
{
	//unregister from GPU
	if (d_matrix != NULL) {
		checkCudaErrors(cudaHostUnregister(this->GetPointerToVoxels()));
		cudaDeviceSynchronize();
	}

	if (this->_matrix != NULL) {
		Deallocate<VoxelType>(this->_matrix);
		this->_matrix = NULL;
	}
	this->_attr._x = 0;
	this->_attr._y = 0;
	this->_attr._z = 0;
	this->_attr._t = 0;
}

template <class VoxelType> VoxelType*  irtkCUGenericImage<VoxelType>::getDevicePoiner() const
{
	return d_matrix;
}


template <class VoxelType> void irtkCUGenericImage<VoxelType>::normalize(VoxelType minVal, VoxelType maxVal) const
{
	normalize_image_gpu(this->GetNumberOfVoxels(), this->getDevicePoiner(), this->getDevicePoiner(), minVal, maxVal);
	checkCudaErrors(cudaDeviceSynchronize());
}

//Using NPP for maximum performance in basic functions
// need to implement each type for NPP
//TODO get min max index and world positions

template <class VoxelType> VoxelType irtkCUGenericImage<VoxelType>::getAverage(int toggle) const
{
#ifdef DUSE_GPL_CODE_ON
	//TODO NPP
	//#ifdef USE_CUDA_NPP
	//could be done with box filter
	//#else
	Reduction<VoxelType> reduction(this->GetPointerToVoxels(), this->GetNumberOfVoxels());
	VoxelType average;
	reduction.Average(&average);
	return average;
	//#endif
#else
	std::cerr << "GPL option not available. Average not computed. NPP version in progress." << std::endl;
	return 0;
#endif
}

#ifdef USE_CUDA_NPP
template <class VoxelType> void irtkCUGenericImage<VoxelType>::GetHistogram(int* hist, int numBins) const
{
	Npp32s *histDevice = nppsMalloc_32s(numBins);
	NppiSize oSizeROI = {this->GetX()*this->GetY(), this->GetZ()*this->GetT()}; // full image
	const int levelCount = numBins + 1;
	VoxelType max_, min_;
	GetMinMax(&min_, &max_);

	computeHistogram_internal(hist, numBins, oSizeROI, histDevice, levelCount, max_, min_);

	nppsFree(histDevice);
	cudaDeviceSynchronize();
}

template<> void irtkCUGenericImage<char>::computeHistogram_internal(int*  hist, int numBins, NppiSize oSizeROI, Npp32s *histDevice, 
	const int levelCount, char max_, char min_) const
{
	printf("ERROR: Histogram for signed char not implemented! No Histogram computed!\n");
}


template<> void irtkCUGenericImage<unsigned char>::computeHistogram_internal(int*  hist, int numBins, NppiSize oSizeROI, Npp32s *histDevice, 
	const int levelCount, unsigned char max_, unsigned char min_) const
{
	int _max = (int)max_;
	int _min = (int)min_;
	int nDeviceBufferSize;
	nppiHistogramRangeGetBufferSize_8u_C1R(oSizeROI, levelCount, &nDeviceBufferSize);
	Npp8u *pDeviceBuffer = nppsMalloc_8u(nDeviceBufferSize);

	Npp32s* levelsHost = new Npp32s[levelCount];
	NPP_CHECK_NPP(nppiEvenLevelsHost_32s(levelsHost, levelCount, min_, max_));

	Npp32s* pLevels = nppsMalloc_32s(levelCount);
	NPP_CHECK_CUDA(cudaMemcpy(pLevels, levelsHost, levelCount * sizeof(Npp32s), cudaMemcpyHostToDevice));

	NPP_CHECK_NPP(nppiHistogramRange_8u_C1R(this->getDevicePoiner(), this->GetX()*this->GetY()*sizeof(Npp8u), oSizeROI,
		histDevice, pLevels, 
		levelCount,
		pDeviceBuffer));

	NPP_CHECK_CUDA(cudaMemcpy(hist, histDevice, numBins * sizeof(Npp32s), cudaMemcpyDeviceToHost));

	//nppsFree(levelsHost);
	nppsFree(pLevels);
	nppsFree(pDeviceBuffer);
}

template<> void irtkCUGenericImage<short>::computeHistogram_internal(int*  hist, int numBins, NppiSize oSizeROI, Npp32s *histDevice, 
	const int levelCount, short max_, short min_) const
{
	int _max = (int)max_;
	int _min = (int)min_;
	int nDeviceBufferSize;
	nppiHistogramRangeGetBufferSize_16s_C1R(oSizeROI, levelCount, &nDeviceBufferSize);
	Npp8u *pDeviceBuffer = nppsMalloc_8u(nDeviceBufferSize);

	Npp32s* levelsHost = new Npp32s[levelCount];
	NPP_CHECK_NPP(nppiEvenLevelsHost_32s(levelsHost, levelCount, min_, max_));

	Npp32s* pLevels = nppsMalloc_32s(levelCount);
	NPP_CHECK_CUDA(cudaMemcpy(pLevels, levelsHost, levelCount * sizeof(Npp32s), cudaMemcpyHostToDevice));

	NPP_CHECK_NPP(nppiHistogramRange_16s_C1R(this->getDevicePoiner(), this->GetX()*this->GetY()*sizeof(Npp16s), oSizeROI,
		histDevice, pLevels, 
		levelCount,
		pDeviceBuffer));

	cudaDeviceSynchronize();
	NPP_CHECK_CUDA(cudaMemcpy(hist, histDevice, numBins * sizeof(Npp32s), cudaMemcpyDeviceToHost));

	nppsFree(pLevels);
	nppsFree(pDeviceBuffer);
	//nppsFree(levelsHost);
}

template<> void irtkCUGenericImage<unsigned short>::computeHistogram_internal(int*  hist, int numBins, NppiSize oSizeROI, Npp32s *histDevice, 
	const int levelCount, unsigned short max_, unsigned short min_) const
{
	int _max = (int)max_;
	int _min = (int)min_;
	int nDeviceBufferSize;
	nppiHistogramRangeGetBufferSize_16u_C1R(oSizeROI, levelCount, &nDeviceBufferSize);
	Npp8u *pDeviceBuffer = nppsMalloc_8u(nDeviceBufferSize);

	Npp32s* levelsHost = new Npp32s[levelCount];
	NPP_CHECK_NPP(nppiEvenLevelsHost_32s(levelsHost, levelCount, min_, max_));

	Npp32s* pLevels = nppsMalloc_32s(levelCount);
	NPP_CHECK_CUDA(cudaMemcpy(pLevels, levelsHost, levelCount * sizeof(Npp32s), cudaMemcpyHostToDevice));

	NPP_CHECK_NPP(nppiHistogramRange_16u_C1R(this->getDevicePoiner(), this->GetX()*this->GetY()*sizeof(Npp16u), oSizeROI,
		histDevice, pLevels, 
		levelCount,
		pDeviceBuffer));

	cudaDeviceSynchronize();
	NPP_CHECK_CUDA(cudaMemcpy(hist, histDevice, numBins * sizeof(Npp32s), cudaMemcpyDeviceToHost));

	//nppsFree(levelsHost);
	nppsFree(pLevels);
	nppsFree(pDeviceBuffer);
}

template<> void irtkCUGenericImage<int>::computeHistogram_internal(int*  hist, int numBins, NppiSize oSizeROI, Npp32s *histDevice, 
	const int levelCount, int max_, int min_) const
{
	printf("sorry, no NPP implementation for 'int' available yet. No histogram computed! \n You could use float.\n");
}

template<> void irtkCUGenericImage<unsigned int>::computeHistogram_internal(int*  hist, int numBins, NppiSize oSizeROI, Npp32s *histDevice, 
	const int levelCount, unsigned int max_, unsigned int min_) const
{
	printf("sorry, no NPP implementation for'unsigned int' available yet. No histogram computed! \n You could use float.\n");
}

template<> void irtkCUGenericImage<float>::computeHistogram_internal(int*  hist, int numBins, NppiSize oSizeROI, Npp32s *histDevice, 
	const int levelCount, float max_, float min_) const
{
	int _max = numeric_limits<int>::max();
	int _min = numeric_limits<int>::min();
	int nDeviceBufferSize;
	nppiHistogramRangeGetBufferSize_32f_C1R(oSizeROI, levelCount, &nDeviceBufferSize);
	Npp8u *pDeviceBuffer = nppsMalloc_8u(nDeviceBufferSize);

	Npp32s* levelsHost = new Npp32s[levelCount];
	NPP_CHECK_NPP(nppiEvenLevelsHost_32s(levelsHost, levelCount, min_, max_));
	//convert to float
	Npp32f* flevelsHost = new Npp32f[levelCount];
	for(int i = 0; i<levelCount; i++){flevelsHost[i] = (float)levelsHost[i]/(float)numeric_limits<int>::max(); printf("%f ",flevelsHost[i]); }

	Npp32f* pLevels = nppsMalloc_32f(levelCount);
	NPP_CHECK_CUDA(cudaMemcpy(pLevels, flevelsHost, levelCount * sizeof(Npp32s), cudaMemcpyHostToDevice));

	NPP_CHECK_NPP(nppiHistogramRange_32f_C1R(this->getDevicePoiner(), this->GetX()*this->GetY()*sizeof(Npp32f), oSizeROI,
		histDevice, pLevels, 
		levelCount,
		pDeviceBuffer));

	cudaDeviceSynchronize();
	NPP_CHECK_CUDA(cudaMemcpy(hist, histDevice, numBins * sizeof(Npp32s), cudaMemcpyDeviceToHost));

	//nppsFree(levelsHost);
	nppsFree(pLevels);
	nppsFree(pDeviceBuffer);
}


template<> void irtkCUGenericImage<double>::computeHistogram_internal(int*  hist, int numBins, NppiSize oSizeROI, Npp32s *histDevice, 
	const int levelCount, double max_, double min_) const
{
	printf("sorry, no NPP implementation for 'double' available yet. No histogram computed! \n You could use float.\n");
}


#endif

template <class VoxelType> void irtkCUGenericImage<VoxelType>::GetMinMax(VoxelType *min, VoxelType *max) const
{
#ifdef USE_CUDA_NPP
	//assuming worst case
	checkCudaErrors(cudaHostRegister(min,sizeof(float), cudaHostRegisterMapped));
	checkCudaErrors(cudaHostRegister(max,sizeof(float), cudaHostRegisterMapped));
	Npp32f* d_min;
	Npp32f* d_max;
	checkCudaErrors(cudaHostGetDevicePointer((void **)&d_min, (void *)min, 0));
	checkCudaErrors(cudaHostGetDevicePointer((void **)&d_max, (void *)max, 0));
	cudaDeviceSynchronize();

	NppiSize oSizeROI = {this->GetX(), this->GetY()*this->GetZ()*this->GetT()};

	// create device scratch buffer for nppiMinMax
	int nDeviceBufferSize;
	NPP_CHECK_NPP(nppiMinMaxGetBufferHostSize_32f_C1R(oSizeROI, &nDeviceBufferSize));

	Npp8u *pDeviceBuffer;
	checkCudaErrors(cudaMalloc((void **)&pDeviceBuffer, nDeviceBufferSize));

	NPP_CHECK_NPP(nppiMinMax_32f_C1R ((Npp32f*)(this->getDevicePoiner()), this->GetX()*sizeof(float), oSizeROI,
		d_min, d_max, pDeviceBuffer));
	cudaDeviceSynchronize();
#else
#ifdef DUSE_GPL_CODE_ON
	Reduction<VoxelType> reduction(this->GetPointerToVoxels(), this->GetNumberOfVoxels());
	reduction.Min(min);
	reduction.Max(max);

#else
	std::cerr << "GPL option not available. NPP not available. No GPU MinMax available! CPU fallback" << std::endl;
	return irtkCUGenericImage<VoxelType>::GetMinMax(min,max);
#endif
#endif
	nppsFree(pDeviceBuffer);
}
#ifdef USE_CUDA_NPP
template <> void irtkCUGenericImage<char>::GetMinMax( char *min,  char *max) const
{
	printf("signed char min max implementation not available yet. using CPU.\n");
	irtkGenericImage<char>::GetMinMax(min, max);
}

template <> void irtkCUGenericImage<short>::GetMinMax( short *min,  short *max) const
{
	Npp16s* d_min = nppsMalloc_16s(1);
	Npp16s* d_max = nppsMalloc_16s(1);

	NppiSize oSizeROI = {this->GetX()*this->GetY(), this->GetZ()*this->GetT()};

	// create device scratch buffer for nppiMinMax
	int nDeviceBufferSize;
	NPP_CHECK_NPP(nppiMinMaxGetBufferHostSize_16s_C1R(oSizeROI, &nDeviceBufferSize));

	Npp8u *pDeviceBuffer;
	checkCudaErrors(cudaMalloc((void **)&pDeviceBuffer, nDeviceBufferSize));

	NPP_CHECK_NPP(nppiMinMax_16s_C1R(this->getDevicePoiner(), this->GetX()*this->GetY()*sizeof(short), oSizeROI,
		d_min, d_max, pDeviceBuffer));

	cudaDeviceSynchronize();
	NPP_CHECK_CUDA(cudaMemcpy(min, d_min, 1 * sizeof(Npp16s), cudaMemcpyDeviceToHost));
	NPP_CHECK_CUDA(cudaMemcpy(max, d_max, 1 * sizeof(Npp16s), cudaMemcpyDeviceToHost));

	nppsFree(pDeviceBuffer);
	nppsFree(d_min);
	nppsFree(d_max);
}

template <> void irtkCUGenericImage<unsigned short>::GetMinMax(unsigned short *min, unsigned short *max) const
{
	checkCudaErrors(cudaHostRegister(min,sizeof(unsigned short), cudaHostRegisterMapped));
	checkCudaErrors(cudaHostRegister(max,sizeof(unsigned short), cudaHostRegisterMapped));
	Npp16u* d_min;
	Npp16u* d_max;
	checkCudaErrors(cudaHostGetDevicePointer((void **)&d_min, (void *)min, 0));
	checkCudaErrors(cudaHostGetDevicePointer((void **)&d_max, (void *)max, 0));
	cudaDeviceSynchronize();

	NppiSize oSizeROI = {this->GetX(), this->GetY()*this->GetZ()*this->GetT()};

	// create device scratch buffer for nppiMinMax
	int nDeviceBufferSize;
	NPP_CHECK_NPP(nppiMinMaxGetBufferHostSize_16u_C1R(oSizeROI, &nDeviceBufferSize));

	Npp8u *pDeviceBuffer;
	checkCudaErrors(cudaMalloc((void **)&pDeviceBuffer, nDeviceBufferSize));

	NPP_CHECK_NPP(nppiMinMax_16u_C1R ((Npp16u*)(this->getDevicePoiner()), this->GetX()*sizeof(unsigned short), oSizeROI,
		d_min, d_max, pDeviceBuffer));
	cudaDeviceSynchronize();

	nppsFree(pDeviceBuffer);
}

template <> void irtkCUGenericImage<unsigned char>::GetMinMax(unsigned char *min, unsigned char *max) const
{
	checkCudaErrors(cudaHostRegister(min,sizeof(unsigned char), cudaHostRegisterMapped));
	checkCudaErrors(cudaHostRegister(max,sizeof(unsigned char), cudaHostRegisterMapped));
	Npp8u* d_min;
	Npp8u* d_max;
	checkCudaErrors(cudaHostGetDevicePointer((void **)&d_min, (void *)min, 0));
	checkCudaErrors(cudaHostGetDevicePointer((void **)&d_max, (void *)max, 0));
	cudaDeviceSynchronize();

	NppiSize oSizeROI = {this->GetX(), this->GetY()*this->GetZ()*this->GetT()};

	// create device scratch buffer for nppiMinMax
	int nDeviceBufferSize;
	NPP_CHECK_NPP(nppiMinMaxGetBufferHostSize_8u_C1R(oSizeROI, &nDeviceBufferSize));

	Npp8u *pDeviceBuffer;
	checkCudaErrors(cudaMalloc((void **)&pDeviceBuffer, nDeviceBufferSize));

	NPP_CHECK_NPP(nppiMinMax_8u_C1R ((Npp8u*)(this->getDevicePoiner()), this->GetX()*sizeof(unsigned char), oSizeROI,
		d_min, d_max, pDeviceBuffer));
	cudaDeviceSynchronize();
	nppsFree(pDeviceBuffer);
}
#endif



template <class VoxelType> void irtkCUGenericImage<VoxelType>::cudaCheckSets()
{
	if(cudaInitDone)
		return;

	cudaDeviceProp deviceProp;
	deviceProp.major = 1;
	deviceProp.minor = 0;
	int minimumComputeVersion = 20;
	int dev = 0; //TODO select device 

	checkCudaErrors(cudaGetDeviceProperties(&deviceProp, dev));

	if ((deviceProp.major * 10 + deviceProp.minor) >= minimumComputeVersion)
	{
		//printf("Using Device %d: %s\n\n", dev, deviceProp.name);
		checkCudaErrors(cudaSetDevice(dev));
	}
	else
	{
		printf("Error: the selected device does not support the minimum compute capability of %d.%d.\n\n",
			minimumComputeVersion / 10, minimumComputeVersion % 10);

		cudaDeviceReset();
		exit(EXIT_FAILURE);
	}

	if (!deviceProp.canMapHostMemory) 
	{
		printf("GPU device does not support page locked mapped memory, which is required for this extension. exit. \n");
		exit(0);
	}
	checkCudaErrors(cudaSetDeviceFlags(cudaDeviceMapHost));

#ifdef USE_CUDA_NPP
	const NppLibraryVersion *libVer   = nppGetLibVersion();
	NppGpuComputeCapability computeCap = nppGetGpuComputeCapability();

	printf("NPP Library Version %d.%d.%d\n", libVer->major, libVer->minor, libVer->build);
#endif

	cudaInitDone = true;
}




/////////implemetations from Template.h to here--

template class irtkCULib_DLLAPI irtkCUGenericImage<char>;
template class irtkCULib_DLLAPI irtkCUGenericImage<unsigned char>;
template class irtkCULib_DLLAPI irtkCUGenericImage<short>;
template class irtkCULib_DLLAPI irtkCUGenericImage<unsigned short>;
template class irtkCULib_DLLAPI irtkCUGenericImage<int>;
template class irtkCULib_DLLAPI irtkCUGenericImage<unsigned int>;
template class irtkCULib_DLLAPI irtkCUGenericImage<float>;
template class irtkCULib_DLLAPI irtkCUGenericImage<double>;

template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkCUGenericImage<unsigned char> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkCUGenericImage<short> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkCUGenericImage<unsigned short> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkCUGenericImage<int> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkCUGenericImage<float> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkCUGenericImage<double> &);

template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkCUGenericImage<char> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkCUGenericImage<short> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkCUGenericImage<unsigned short> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkCUGenericImage<int> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkCUGenericImage<float> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkCUGenericImage<double> &);

template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkCUGenericImage<char> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkCUGenericImage<unsigned char> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkCUGenericImage<unsigned short> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkCUGenericImage<int> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkCUGenericImage<float> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkCUGenericImage<double> &);

template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkCUGenericImage<char> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkCUGenericImage<unsigned char> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkCUGenericImage<short> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkCUGenericImage<int> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkCUGenericImage<float> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkCUGenericImage<double> &);

template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkCUGenericImage<char> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkCUGenericImage<unsigned char> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkCUGenericImage<short> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkCUGenericImage<unsigned short> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkCUGenericImage<float> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkCUGenericImage<double> &);

template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkCUGenericImage<char> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkCUGenericImage<unsigned char> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkCUGenericImage<short> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkCUGenericImage<unsigned short> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkCUGenericImage<int> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkCUGenericImage<double> &);

template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkCUGenericImage<char> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkCUGenericImage<unsigned char> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkCUGenericImage<short> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkCUGenericImage<unsigned short> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkCUGenericImage<int> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkCUGenericImage<float> &);

//template irtkCUGenericImage<float>& irtkCUGenericImage<float>::operator=(const irtkCUGenericImage<short> &);


template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkGenericImage<unsigned char> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkGenericImage<short> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkGenericImage<unsigned short> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkGenericImage<int> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkGenericImage<float> &);
template irtkCUGenericImage<char>::irtkCUGenericImage(const irtkGenericImage<double> &);

template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkGenericImage<char> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkGenericImage<short> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkGenericImage<unsigned short> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkGenericImage<int> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkGenericImage<float> &);
template irtkCUGenericImage<unsigned char>::irtkCUGenericImage(const irtkGenericImage<double> &);

template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkGenericImage<char> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkGenericImage<unsigned char> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkGenericImage<unsigned short> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkGenericImage<int> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkGenericImage<float> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkGenericImage<double> &);
template irtkCUGenericImage<short>::irtkCUGenericImage(const irtkGenericImage<short> &);

template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkGenericImage<char> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkGenericImage<unsigned char> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkGenericImage<short> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkGenericImage<int> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkGenericImage<float> &);
template irtkCUGenericImage<unsigned short>::irtkCUGenericImage(const irtkGenericImage<double> &);

template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkGenericImage<char> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkGenericImage<unsigned char> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkGenericImage<short> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkGenericImage<unsigned short> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkGenericImage<float> &);
template irtkCUGenericImage<int>::irtkCUGenericImage(const irtkGenericImage<double> &);

template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkGenericImage<char> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkGenericImage<unsigned char> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkGenericImage<short> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkGenericImage<unsigned short> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkGenericImage<int> &);
template irtkCUGenericImage<float>::irtkCUGenericImage(const irtkGenericImage<double> &);

template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkGenericImage<char> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkGenericImage<unsigned char> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkGenericImage<short> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkGenericImage<unsigned short> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkGenericImage<int> &);
template irtkCUGenericImage<double>::irtkCUGenericImage(const irtkGenericImage<float> &);
