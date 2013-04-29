Project(irtkCUcommonLib)

cmake_minimum_required(VERSION 2.6)

OPTION(BUILD_WITH_GPL_CODE "Use also GPL fallback code. Do not use for public releases!!" OFF)
if(BUILD_WITH_GPL_CODE)
 add_definitions(-DUSE_GPL_CODE_ON)
endif(BUILD_WITH_GPL_CODE)

#if IRTK_SOURCE_DIR is set, we are building within the IRTK source tree
if(IRTK_SOURCE_DIR)
SET(BUILD_OUT_OF_IRTK_TREE FALSE)
else(IRTK_SOURCE_DIR)
SET(BUILD_OUT_OF_IRTK_TREE TRUE)
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

# CUDA
option(CUDA_BUILD_EMULATION "enable emulation mode" OFF)
find_package(CUDA REQUIRED)
# FindCUDA.cmake of CMake 2.8.7 does not look for NPP
if (NOT DEFINED CUDA_npp_LIBRARY AND COMMAND find_cuda_helper_libs AND NOT CUDA_VERSION VERSION_LESS "4.0")
  find_cuda_helper_libs(npp)
endif ()
include_directories(${CUDA_INCLUDE_DIRS})
# hide variable as SDK now (usually) part of CUDA Toolkit
mark_as_advanced(FORCE CUDA_SDK_ROOT_DIR)
# additional directories where to look for SDK which are
# not yet added in FindCUDA.cmake of CMake
list(APPEND CUDA_SDK_SEARCH_PATH "${CUDA_TOOLKIT_ROOT_DIR}/samples") 
# look for required SDK header files as suggested in FindCUDA.cmake
find_path(CUDA_HELPER_INCLUDE_DIR
  helper_cuda.h
  PATHS         ${CUDA_SDK_SEARCH_PATH}
  PATH_SUFFIXES "common/inc"
  DOC           "Location of helper_cuda.h of the CUDA SDK."
  NO_DEFAULT_PATH
)
find_path(CUDA_HELPER_INCLUDE_DIR helper_cuda.h DOC "Location of helper_cuda.h of the CUDA SDK.")
mark_as_advanced(FORCE CUDA_HELPER_INCLUDE_DIR)
# either add found path to include search path or show user the
# appropriate CMake variables to tell us where to find the headers
if (CUDA_HELPER_INCLUDE_DIR)
  include_directories(${CUDA_HELPER_INCLUDE_DIR})
elseif (NOT CUDA_SDK_ROOT_DIR)
  mark_as_advanced(CLEAR CUDA_SDK_ROOT_DIR)
else ()
  mark_as_advanced(CLEAR CUDA_HELPER_INCLUDE_DIR)
endif ()

if( ${CUDA_VERSION} VERSION_LESS "5")
  message(INFO "Building with CUDA < 5.0")
else()
  message(STATUS "Building with CUDA >= 5.0")
  # we build for all useful compute capabilities (C.P > 2.0)
  set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-gencode;arch=compute_12,code=sm_12;-gencode;arch=compute_20,code=sm_21;-gencode;arch=compute_20,code=sm_20;-gencode;arch=compute_30,code=sm_30;-use_fast_math;")
  #untested for CP < 2.0
  #set(CUDA_NVCC_FLAGS "-gencode;arch=compute_12,code=sm_12;-use_fast_math;")
endif()

# check if we have c++11
set(CMAKE_CXX_COMPILER_TMP ${CMAKE_CXX_COMPILER})
set(CMAKE_CXX_COMPILER ${CUDA_NVCC_EXECUTABLE})
include (cmake/CheckCxxFeatures.cmake)
cxx_check_feature ("c++11" "nullptr" "N2431" HAVE_CXX11_NULLPTR "${ADDITIONAL_DEFINITIONS}")
set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER_TMP})

if(HAVE_CXX11_NULLPTR)
  add_definitions(-DHAVE_CXX11_NULLPTR)
else(HAVE_CXX11_NULLPTR)
  message(WARNING "Your nvcc does not support c++11 yet. Using own nullptr implementation!")
endif(HAVE_CXX11_NULLPTR)


if(BUILD_OUT_OF_IRTK_TREE)

# we do not need to search for IRTK for in source tree build
find_package(IRTK REQUIRED)
message("${IRTK_LIBRARIES_DIR}")
include_directories(${IRTK_INCLUDE_DIRS})

else(BUILD_OUT_OF_IRTK_TREE)

include_directories(../common++/include ../image++/include ../geometry++/include ../contrib++/include)

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
)

if(BUILD_WITH_GPL_CODE)
set(irtkCUcommonLib_SRC ${irtkCUcommonLib_SRC} 
			src/reduction/MaxKernel.cu
			src/reduction/MinKernel.cu
			src/reduction/reduction.cpp
			src/reduction/SumKernel.cu
)
endif(BUILD_WITH_GPL_CODE)

set(irtkCUcommonLib_HDRS include/irtkCUGenericImage.h
			 include/irtkCUImage.h
			 src/IO/irtkCUImage_kernels.cuh
			 include/irtkCUSharedLibMacros.h
			 include/irtkCUAbstractFilterBase.h
			 include/irtkCUMaskFilter.h
)

if(BUILD_WITH_GPL_CODE)
set(irtkCUcommonLib_HDRS ${irtkCUcommonLib_HDRS} 
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
endif(BUILD_WITH_GPL_CODE)

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



