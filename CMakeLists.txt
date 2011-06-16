cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  # policy settings ...
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

PROJECT(IRTK)

INCLUDE(CTest)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Wno-deprecated")        ## Optimize
endif()

OPTION(USE_VXL "Use VXL library instead of Numerical Recipes" OFF)
IF (USE_VXL)
  SUBDIRS(vxl)
  #  INCLUDE_DIRECTORIES(${VCL_INCLUDE_DIR} ${VXLCORE_INCLUDE_DIR})
  INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/vxl/core ${IRTK_SOURCE_DIR}/vxl/vcl)
  ADD_DEFINITIONS(-DUSE_VXL)
ELSE (USE_VXL)
  SUBDIRS(recipes/src)
ENDIF (USE_VXL)

SUBDIRS(common++/src geometry++/src image++/src contrib++/src applications packages documentation)

INCLUDE(${IRTK_SOURCE_DIR}/config/IRTKInclude.cmake)

SET(EXECUTABLE_OUTPUT_PATH ${IRTK_BINARY_DIR}/bin CACHE STRING
  "Where to put the executables")
SET(LIBRARY_OUTPUT_PATH ${IRTK_BINARY_DIR}/lib CACHE STRING
  "Where to put the libraries")

IF (BUILD_WITH_NIFTI)
   SUBDIRS(nifti)
ENDIF (BUILD_WITH_NIFTI)

IF (ENABLE_WRAPPING)
  SUBDIRS(wrapping)
ENDIF (ENABLE_WRAPPING)

CONFIGURE_FILE(${IRTK_SOURCE_DIR}/config/IRTKInclude.cmake
               ${IRTK_BINARY_DIR}/lib/IRTKInclude.cmake)	

CONFIGURE_FILE(${IRTK_SOURCE_DIR}/InstallerScripts/irtkInstaller.sh
               ${IRTK_BINARY_DIR}/bin/irtkInstaller.sh)
               
EXEC_PROGRAM(chmod ${IRTK_BINARY_DIR}/bin ARGS +x;irtkInstaller.sh)
