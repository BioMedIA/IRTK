cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  # policy settings ...
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

PROJECT(IRTK)

INCLUDE(CTest)

OPTION(USE_VXL "Use VXL library instead of Numerical Recipes" OFF)
IF (USE_VXL)
  SUBDIRS(vxl)
  #  INCLUDE_DIRECTORIES(${VCL_INCLUDE_DIR} ${VXLCORE_INCLUDE_DIR})
  INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/vxl/core ${IRTK_SOURCE_DIR}/vxl/vcl)
  ADD_DEFINITIONS(-DUSE_VXL)
ELSE (USE_VXL)
  SUBDIRS(recipes/src)
ENDIF (USE_VXL)

OPTION(BUILD_WITH_TBB "Use TBB library for multithreading support" OFF)

OPTION(BUILD_GPU_SUPPORT_WITH_CUDA "Use a preliminary GPU accelerated version of IRTK" OFF)
IF (BUILD_GPU_SUPPORT_WITH_CUDA)
  add_subdirectory(cudaCommon++)
  add_subdirectory(packages/segmentation/src/irtkReconstructionCuda)
ENDIF(BUILD_GPU_SUPPORT_WITH_CUDA)

SUBDIRS(common++/src geometry++/src image++/src contrib++/src applications packages documentation wrapping)

INCLUDE(${IRTK_SOURCE_DIR}/config/IRTKInclude.cmake)

SET(EXECUTABLE_OUTPUT_PATH ${IRTK_BINARY_DIR}/bin CACHE STRING
  "Where to put the executables")
SET(LIBRARY_OUTPUT_PATH ${IRTK_BINARY_DIR}/lib CACHE STRING
  "Where to put the libraries")

IF (BUILD_WITH_NIFTI)
   SUBDIRS(nifti)
ENDIF (BUILD_WITH_NIFTI)

CONFIGURE_FILE(${IRTK_SOURCE_DIR}/config/IRTKInclude.cmake
               ${IRTK_BINARY_DIR}/lib/IRTKInclude.cmake)	

CONFIGURE_FILE(${IRTK_SOURCE_DIR}/InstallerScripts/irtkInstaller.sh
               ${IRTK_BINARY_DIR}/bin/irtkInstaller.sh)
               
EXEC_PROGRAM(chmod ${IRTK_BINARY_DIR}/bin ARGS +x;irtkInstaller.sh)

INCLUDE(${IRTK_SOURCE_DIR}/config/IRTKGenerateConfigFile.cmake)
