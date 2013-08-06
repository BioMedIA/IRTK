cmake_minimum_required(VERSION 2.6)

# CMakeLists.txt adapted from
# https://github.com/lantiga/cyitk

#project(CYTHON)

find_package(PythonLibs REQUIRED)
if(USE_PYTHONLIBS_FILE)
  include(${USE_PYTHONLIBS_FILE})
endif(USE_PYTHONLIBS_FILE)

find_package(PythonInterp REQUIRED)
if(USE_PYTHONINTERP_FILE)
  include(${USE_PYTHONINTERP_FILE})
endif(USE_PYTHONINTERP_FILE)

#get numpy include path
execute_process(
  COMMAND
  ${PYTHON_EXECUTABLE} -c "import numpy; print numpy.get_include()"
  OUTPUT_VARIABLE NUMPY_INCLUDE_PATH
  RESULT_VARIABLE NUMPY_ERR
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
if(NUMPY_ERR)
  message(SEND_ERROR "WARNING: NumPy header not found.")
endif(NUMPY_ERR)
set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})

#TODO: Windows? Need a FindCython.cmake
find_program(CYTHON_EXECUTABLE cython /usr/bin /usr/local/bin $ENV{PYTHONPATH} ~/.local/bin)

#create temporary directory
set(TMP_SRC_DIR ${CMAKE_BINARY_DIR}/tmp)
file(MAKE_DIRECTORY ${TMP_SRC_DIR})

#generate templates
execute_process(
  COMMAND
  ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/wrapping/cython/generate_templates.py ${TMP_SRC_DIR} ${CMAKE_SOURCE_DIR}/wrapping/cython/
  OUTPUT_VARIABLE PYTHON_OUTPUT
  ERROR_VARIABLE PYTHON_ERROR
  RESULT_VARIABLE PYTHON_ERR
  )
if(PYTHON_OUTPUT)
  message(${PYTHON_OUTPUT})
endif(PYTHON_OUTPUT)
if(PYTHON_ERR)
  message(SEND_ERROR ${PYTHON_ERROR})
endif(PYTHON_ERR)

set( CPP_SRCS
  "src/registration.cc"
  "src/reconstruction.cc"
  "${TMP_SRC_DIR}/templates.cc"
  "src/irtk2cython.cc"
  "src/voxellise.cc"
  "src/drawing.cc"
  )

include_directories(${PYTHON_INCLUDE_PATH})
include_directories("include")
include_directories(${CMAKE_BINARY_DIR}"/tmp")

#run cython on all pyx files
execute_process(
  COMMAND
  ${CYTHON_EXECUTABLE} --cplus ${CMAKE_BINARY_DIR}/tmp/_irtk.pyx -o ${CMAKE_BINARY_DIR}/tmp/_irtk.cpp
  OUTPUT_VARIABLE CYTHON_OUTPUT
  ERROR_VARIABLE CYTHON_ERROR
  RESULT_VARIABLE CYTHON_ERR
  )
if(CYTHON_ERR)
  message(SEND_ERROR ${CYTHON_ERROR})
endif(CYTHON_ERR)
  
set(CPP_SRCS ${CPP_SRCS} ${CMAKE_BINARY_DIR}/tmp/_irtk.cpp)
  
add_library( "_irtk" MODULE ${CPP_SRCS})
target_link_libraries( "_irtk" ${IRTK_LIBRARIES} ${PYTHON_LIBRARY})
set_target_properties("_irtk" PROPERTIES PREFIX "" LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/irtk)

# copy all python files
file( COPY ${CMAKE_SOURCE_DIR}/wrapping/cython/irtk DESTINATION ${CMAKE_BINARY_DIR}/lib )

