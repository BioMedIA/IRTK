# Compatible with CMake 2.6.x
GET_FILENAME_COMPONENT(CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)

# Prefer project-own Find<Package>.cmake modules for FIND_PACKAGE().
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

# Option to produce condor executables
OPTION(BUILD_CONDOR_EXE "Build condor executables" OFF)
MARK_AS_ADVANCED(BUILD_CONDOR_EXE)
IF (BUILD_CONDOR_EXE)
   SET (CMAKE_X_LIBS "-lGL -L/usr/X11R6/lib -lm -lXext -lXt -lSM -lX11  -lICE -ldl -lnsl")
   SET (CMAKE_MODULE_LINK_FLAGS "-static")
   SET (CMAKE_SHLIB_LINK_FLAGS "-static")
ENDIF (BUILD_CONDOR_EXE)

# Option to build with OpenCV or not.
OPTION(BUILD_WITH_OPENCV "Build using OPENCV" OFF)

IF (BUILD_WITH_OPENCV)
  FIND_PACKAGE(OpenCV REQUIRED)
  ADD_DEFINITIONS(-DHAS_OPENCV)
  INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIR})
  LINK_LIBRARIES(${OpenCV_LIBS})
ENDIF (BUILD_WITH_OPENCV)

# Option to produce multi-threaded executables using TBB
OPTION(BUILD_WITH_TBB "Build multi-threaded executables using TBB" OFF)

IF (BUILD_WITH_TBB)
  FIND_PACKAGE(TBB REQUIRED)
  IF (TBB_FOUND)
    # Attention: DO NOT define TBB_DEPRECATED by default or before including the
    #            other TBB header files, in particular parallel_for. The deprecated
    #            behavior of parallel_for is to not choose the chunk size (grainsize)
    #            automatically!
    #
    # http://software.intel.com/sites/products/documentation/doclib/tbb_sa/help/tbb_userguide/Automatic_Chunking.htm
    ADD_DEFINITIONS(-DHAS_TBB)
    INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIRS})
    LINK_DIRECTORIES(${TBB_LIBRARY_DIRS})
    LINK_LIBRARIES(${TBB_LIBRARIES})
  ENDIF (TBB_FOUND)
ENDIF(BUILD_WITH_TBB)

INCLUDE(${CMAKE_ROOT}/Modules/FindZLIB.cmake)

IF (ZLIB_FOUND)
  ADD_DEFINITIONS(-DHAS_ZLIB -DHAVE_ZLIB)
  INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
  LINK_LIBRARIES(${ZLIB_LIBRARIES})
ENDIF (ZLIB_FOUND)

INCLUDE(${CMAKE_ROOT}/Modules/FindPNG.cmake)

FIND_PACKAGE(FLTK REQUIRED)

IF (FLTK_FOUND)
  INCLUDE_DIRECTORIES(${FLTK_INCLUDE_DIR})
  LINK_LIBRARIES (${FLTK_LIBRARIES})
ENDIF (FLTK_FOUND)

IF (APPLE)
  FIND_PROGRAM(APPLE_REZ Rez /Developer/Tools)
  FIND_FILE(FLTK_REZ_FILE mac.r ${FLTK_INCLUDE_DIR})
  MARK_AS_ADVANCED(APPLE_REZ)
  MARK_AS_ADVANCED(FLTK_REZ_FILE)
ENDIF (APPLE)

IF (NOT BUILD_CONDOR_EXE)

   INCLUDE (${CMAKE_ROOT}/Modules/FindGLUT.cmake)

   IF (GLUT_INCLUDE_PATH)
      INCLUDE_DIRECTORIES(${GLUT_INCLUDE_PATH})
   ENDIF(GLUT_INCLUDE_PATH)

   IF (GLUT_LIBRARY)
      LINK_LIBRARIES (${GLUT_LIBRARY})
   ENDIF (GLUT_LIBRARY)

   INCLUDE (${CMAKE_ROOT}/Modules/FindOpenGL.cmake)

   IF (GLU_LIBRARY)
      LINK_LIBRARIES (${GLU_LIBRARY})
   ENDIF (GLU_LIBRARY)
 
   IF (OPENGL_INCLUDE_PATH)
      INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_PATH})
   ENDIF(OPENGL_INCLUDE_PATH)

   IF (OPENGL_LIBRARY)
      LINK_LIBRARIES (${OPENGL_LIBRARY})
   ENDIF (OPENGL_LIBRARY)

ENDIF (NOT BUILD_CONDOR_EXE)

INCLUDE (${CMAKE_ROOT}/Modules/Dart.cmake)

IF (WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /TP /Ze /W0")
  ADD_DEFINITIONS(-DvtkCommon_EXPORTS)
ELSE (WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-deprecated -Wno-write-strings")
ENDIF (WIN32)

# Option to wrap python.
OPTION(WRAP_PYTHON "Generate python wrapper libraries" OFF)

IF (WRAP_PYTHON)
  SET(ENABLE_WRAPPING TRUE)
  INCLUDE(${CMAKE_ROOT}/Modules/FindPythonLibs.cmake)
ENDIF (WRAP_PYTHON)

# Find SWIG.
IF (ENABLE_WRAPPING)
  INCLUDE(${CMAKE_ROOT}/Modules/FindSWIG.cmake)
  IF (SWIG_FOUND)
    INCLUDE(${SWIG_USE_FILE})  
    IF (UNIX)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
      SET(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    ENDIF (UNIX)
  ENDIF (SWIG_FOUND)
ENDIF (ENABLE_WRAPPING)

OPTION(WRAP_CYTHON "Generate cython wrapper libraries" OFF)
IF (WRAP_CYTHON)
  IF (UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  ENDIF (UNIX)
ENDIF (WRAP_CYTHON)

ADD_DEFINITIONS(-DIMPERIAL -DANSI -DHAS_CONTRIB -DNO_BOUNDS -DENABLE_UNIX_COMPRESS)

INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/recipes/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/common++/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/geometry++/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/image++/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/contrib++/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/packages/transformation/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/packages/registration/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/packages/registration2/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/packages/segmentation/include)
INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/external/gco-v3.0)

LINK_DIRECTORIES(${IRTK_BINARY_DIR}/lib) 

# Option to build with PNG or not.
OPTION(BUILD_WITH_PNG "Build using PNG" OFF)

IF (BUILD_WITH_PNG)
  ADD_DEFINITIONS(-DHAS_PNG)
  INCLUDE_DIRECTORIES(${PNG_INCLUDE_DIR})
  LINK_LIBRARIES(${PNG_LIBRARIES})
ENDIF (BUILD_WITH_PNG)

# Option to build with VTK or not.
OPTION(BUILD_WITH_VTK "Build using VTK" OFF)

IF (BUILD_WITH_VTK)
   # Add VTK
   INCLUDE(${CMAKE_ROOT}/Modules/FindVTK.cmake)

   IF (VTK_FOUND)
      ADD_DEFINITIONS(-DHAS_VTK)
      INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
      LINK_DIRECTORIES(${VTK_LIBRARY_DIRS})

      # Add patented library if available
      IF (VTK_KITS MATCHES "PATENTED")
         ADD_DEFINITIONS(-DHAS_VTK_PATENTED)
	  LINK_LIBRARIES (vtkPatented)
      ENDIF (VTK_KITS MATCHES "PATENTED")

       # Add patented library if available
      IF (VTK_KITS MATCHES "HYBRID")
         ADD_DEFINITIONS(-DHAS_VTK_HYBRID)
	 LINK_LIBRARIES (vtkHybrid)
      ENDIF (VTK_KITS MATCHES "HYBRID")

     LINK_LIBRARIES (vtkRendering vtkImaging
      vtkGraphics vtkFiltering vtkIO vtkCommon)
   ENDIF (VTK_FOUND)
ENDIF (BUILD_WITH_VTK)


LINK_LIBRARIES(segmentation++ registration2++ registration++ transformation++ contrib++
image++ geometry++ common++ recipes) 

# Options to build with nifti, znz and possibly fslio
OPTION(BUILD_WITH_NIFTI "Build using NIFTI support" ON)
IF (BUILD_WITH_NIFTI)
   ADD_DEFINITIONS(-DHAS_NIFTI)
   INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/nifti/niftilib)
   INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/nifti/znzlib)
   LINK_LIBRARIES(znz)
   LINK_LIBRARIES(niftiio)
ENDIF (BUILD_WITH_NIFTI)

# Option to build with cardiac spatial temporal correction, segmentation and motion tracking toolbox.
OPTION(BUILD_CARDIAC "Build with cardiac tool box" OFF)
IF (BUILD_CARDIAC)
   ADD_DEFINITIONS(-DHAS_CARDIAC)
   INCLUDE_DIRECTORIES(${IRTK_SOURCE_DIR}/packages/cardiac/include)
   LINK_LIBRARIES(cardiac++)
ENDIF (BUILD_CARDIAC)

