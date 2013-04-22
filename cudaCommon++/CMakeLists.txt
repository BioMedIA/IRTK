Project(irtkCUcommonLib)

cmake_minimum_required(VERSION 2.6)

#if IRTK_SOURCE_DIR is set, we are building within the IRTK source tree
#but we leave it as option in case somewone wants to do something strange
if(IRTK_SOURCE_DIR)
SET(BUILD_OUT_OF_IRTK_TREE TRUE)
#OPTION(BUILD_OUT_OF_IRTK_TREE "build outside of the IRTK main source tree." ON)
else(IRTK_SOURCE_DIR)
SET(BUILD_OUT_OF_IRTK_TREE FALSE)
#OPTION(BUILD_OUT_OF_IRTK_TREE "build outside of the IRTK main source tree." OFF)
endif(IRTK_SOURCE_DIR)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ./cmake)

if (WIN32)
  # for msvc 8 deprecated warning removal
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
  add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE_GLOBALS)
  add_definitions(-DWIN32 -D_WIN32)
  add_definitions(-DNOMINMAX)

  #set(CMAKE_DEBUG_POSTFIX "D")

  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /NODEFAULTLIB:LIBCMT.lib;libcmtd.lib")
endif(WIN32)

if (WIN32)
  add_definitions(-DGL_GLEXT_PROTOTYPES -bigobj)
endif()

# Cuda:
set(CUDA_BUILD_EMULATION OFF CACHE BOOL "enable emulation mode")
find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})
#todo make smarter
include_directories(${CUDA_SDK_ROOT_DIR}/common/inc)

if( ${CUDA_VERSION} VERSION_LESS "5")
MESSAGE(INFO "Building with CUDA < 5.0")
message( SEND_ERROR "Softshell requires CUDA 5.0 or higher" )
else()
MESSAGE(INFO "Building with CUDA >= 5.0")
# we build for all useful compute capabilities (C.P > 2.0)
set(CUDA_NVCC_FLAGS "-gencode;arch=compute_20,code=sm_21;-gencode;arch=compute_20,code=sm_20;-gencode;arch=compute_30,code=sm_30;-use_fast_math;")
#untested for CP < 2.0
#set(CUDA_NVCC_FLAGS "-gencode;arch=compute_12,code=sm_12;-use_fast_math;")
endif()

include_directories(${CUDA_INCLUDE_DIRS} ${CUDA_SDK_ROOT_DIR}/common/inc)

if(BUILD_OUT_OF_IRTK_TREE)

# we do not need to search for IRTK for in source tree build
find_package(IRTK REQUIRED)
message("${IRTK_LIBRARIES_DIR}")
include_directories(${IRTK_INCLUDE_DIRS})

else(BUILD_OUT_OF_IRTK_TREE)

include_directories(../common++/include ../image++/include)

SET(IRTK_LIBRARIES contrib++ image++ geometry++ common++) 
IF (USE_VXL)
  SET(IRTK_LIBRARIES irtkvnl irtkvcl irtkvnl_algo irtknetlib ${IRTK_LIBRARIES})
ELSE (USE_VXL)
  SET(IRTK_LIBRARIES recipes ${IRTK_LIBRARIES})
ENDIF (USE_VXL)

endif(BUILD_OUT_OF_IRTK_TREE)


include_directories(include)

set(irtkCUcommonLib_SRC src/IO/irtkCUGenericImage.cpp
			src/IO/irtkCUImage_kernels.cu
			src/filters/irtkCUMaskFilter.cpp
			src/reduction/MaxKernel.cu
			src/reduction/MinKernel.cu
			src/reduction/reduction.cpp
			src/reduction/SumKernel.cu
)
set(irtkCUcommonLib_HDRS include/irtkCUGenericImage.h
			 include/irtkCUImage.h
			 src/IO/irtkCUImage_kernels.cuh
			 include/irtkCUSharedLibMacros.h
			 include/irtkCUAbstractFilterBase.h
			 include/irtkCUMaskFilter.h
			 include/memory/BaseArray.h
			 include/memory/BaseMatrix.h
			 include/memory/CudaArray.h
			 include/memory/CudaMatrix.h
			 include/memory/DeviceAccessibleVariable.h
			 include/memory/DeviceArray.h
			 include/memory/DeviceMatrix.h
			 include/memory/HostArray.h
			 include/memory/HostMatrix.h
			 include/reduction/reduction.h
			 include/reduction/SumWarp.h
			 include/reduction/Utilities.h
			 include/reduction/CudaStreams.h
			 include/reduction/CudaDefinitions.h
)

if(CUDA_npp_LIBRARY)
message("building for NPP library support in ${CUDA_npp_LIBRARY}") 
add_definitions(-DUSE_CUDA_NPP)
endif(CUDA_npp_LIBRARY)

link_directories(${IRTK_LIBRARIES_DIR})

#TODO make IRTK SHARED
add_definitions(-DirtkCUcommonLib_USE_STATIC)
cuda_add_library(irtkCUcommonLib STATIC ${irtkCUcommonLib_SRC} ${irtkCUcommonLib_HDRS})
target_link_libraries(irtkCUcommonLib ${IRTK_LIBRARIES} ${CUDA_npp_LIBRARY})

if(NOT BUILD_OUT_OF_IRTK_TREE)
INSTALL_TARGETS(/lib irtkCUcommonLib)
endif(NOT BUILD_OUT_OF_IRTK_TREE)

#tests
if(BUILD_TESTING)
cuda_add_executable(irtkCUImageTest tests/irtkCUImageTest.cpp)
target_link_libraries(irtkCUImageTest irtkCUcommonLib ${IRTK_LIBRARIES})

if(NOT BUILD_OUT_OF_IRTK_TREE)
INSTALL_TARGETS(/bin irtkCUImageTest)
endif(NOT BUILD_OUT_OF_IRTK_TREE)
endif(BUILD_TESTING)



