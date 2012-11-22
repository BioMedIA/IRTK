# Copied on 19 October 2012 from
# https://oceanai.mit.edu/svn/moos-ivp-aro/releases/moos-ivp-4.2/MOOS/cmake_modules/FindFLTK.cmake
#
# Change history:
#
# 22/10/2012 Andreas Schuh   - Do not attempt to use CMake cache entries from FLTK build tree.
#                            - Do not set FLTK_INCLUDE_DIR when using FLTKConfig.cmake.
#                            - Use proper search paths to find FLTK_DIR.
#                            - Acknowlegde FLTK_DIR when looking for non-CMake FLTK installation
#                              by only looking for files within the specified directory.
#                            - Slight change of verbose output messages.
#                            - Use X11_Xft_LIB variable set by FindX11.cmake module of CMake.
#                            - Fix raising of error when FLTK not found by changing FLTK_REQUIRED to FLTK_FIND_REQUIRED.
#                            - Change FLTK_IMAGES_LIBS to FLTK_IMAGES_LIBRARIES.
#                            - Remove backwards compatible settings for CMake 1.4.
# 14/11/2012 Andreas Schuh   - Set FLTK_INCLUDE_DIR to FLTK_INCLUDE_DIRS as used by FLTK 1.1.
# 22/11/2012 Andreas Schuh   - Fix/improve support of FLTK build and installed on Windows.
#
# Find the native FLTK includes and library.  This is a hacked about version of
# the FindFLTK.cmake module provided with CMake.
#
# The following settings are defined.
# FLTK_FOUND            = whether FLTK was found and the following variables are valid
# FLTK_INCLUDE_DIR      = where to find include files
# FLTK_LIBRARIES        = list of fltk libraries
# FLTK_FLUID_EXECUTABLE = full path to the fluid command
# FLTK_WRAP_UI          = this allows the FLTK_WRAP_UI command to work. True, if fluid executable found

# The following settings should not be used in general.
# FLTK_BASE_LIBRARY    = the full path to fltk library
# FLTK_GL_LIBRARY      = the full path to fltk_gl library
# FLTK_FORMS_LIBRARY   = the full path to fltk_forms library
# FLTK_IMAGES_LIBRARY  = the full path to fltk_images library

# Help string for FLTK_DIR.
SET(_FLTK_DIR_STRING "directory containing the FLTKConfig.cmake file, e.g.,\n\t")
IF(WIN32)
  SET(_FLTK_DIR_STRING "${_FLTK_DIR_STRING}C:\\Program Files\\FLTK\\CMake")
ELSE(WIN32)
  SET(_FLTK_DIR_STRING "${_FLTK_DIR_STRING}/usr/local/lib/FLTK-1.3")
ENDIF(WIN32)
SET(_FLTK_DIR_STRING "${_FLTK_DIR_STRING}\nor installation prefix of FLTK")

# Platform dependent libraries required by FLTK
IF(WIN32)
  IF(NOT CYGWIN)
    IF(BORLAND)
      SET(_FLTK_PLATFORM_DEPENDENT_LIBRARIES import32)
    ELSE(BORLAND)
      SET(_FLTK_PLATFORM_DEPENDENT_LIBRARIES wsock32 comctl32)
    ENDIF(BORLAND)
  ENDIF(NOT CYGWIN)
ENDIF(WIN32)

IF(UNIX)
  FIND_PACKAGE(X11 QUIET)
  IF(X11_FOUND)
    SET(_FLTK_PLATFORM_DEPENDENT_LIBRARIES ${X11_LIBRARIES} -lm)
  ENDIF(X11_FOUND)
ENDIF(UNIX)

IF(APPLE)
  SET(_FLTK_PLATFORM_DEPENDENT_LIBRARIES "-framework Carbon -framework AGL -framework Cocoa -framework ApplicationServices -lz")
ENDIF(APPLE)

IF(CYGWIN)
  SET(_FLTK_PLATFORM_DEPENDENT_LIBRARIES ole32 uuid comctl32 wsock32 supc++ -lm -lgdi32)
ENDIF(CYGWIN)

# Initialize FLTK_DIR if not already known.
IF(NOT FLTK_DIR)
  # Get the system search path for executable files as a CMake list.
  IF(UNIX)
    STRING(REGEX MATCHALL "[^:]+" _FLTK_DIR_SEARCH_PATH "$ENV{PATH}")
  ELSE(UNIX)
    STRING(REGEX REPLACE "\\\\" "/" _FLTK_DIR_SEARCH_PATH "$ENV{PATH}")
  ENDIF(UNIX)
  STRING(REGEX REPLACE "/(;|$)" ";" _FLTK_DIR_SEARCH_PATH "${_FLTK_DIR_SEARCH_PATH}")

  # Construct a set of paths relative to the system search path for executable files.
  SET(_FLTK_CONFIG_SEARCH "")
  SET(_FLTK_DIR_SEARCH    "")
  FOREACH(_FLTK_dir ${_FLTK_DIR_SEARCH_PATH})
    IF(WIN32)
      SET(_FLTK_CONFIG_SEARCH ${_FLTK_CONFIG_SEARCH} "${_FLTK_dir}/.." "${_FLTK_dir}/../CMake")
    ELSE(WIN32)
      SET(_FLTK_CONFIG_SEARCH ${_FLTK_CONFIG_SEARCH} "${_FLTK_dir}/.." "${_FLTK_dir}/../lib/fltk" "${_FLTK_dir}/../lib/FLTK")
    ENDIF(WIN32)
    SET(_FLTK_DIR_SEARCH ${_FLTK_DIR_SEARCH} "${_FLTK_dir}/..")
  ENDFOREACH(_FLTK_dir)
  UNSET(_FLTK_dir)

  # Look for an installation or build tree. 
  FIND_PATH(FLTK_DIR
    NAMES FLTKConfig.cmake
    PATHS
      # Look for an environment variable FLTK_DIR.
      $ENV{FLTK_DIR}
      # Look in places relative to the system executable search path.
      ${_FLTK_CONFIG_SEARCH}
      # Look in standard UNIX install locations.
      /usr/local/lib/fltk
      /usr/lib/fltk
      /usr/local/include
      /usr/include
      /usr/local/fltk
      /usr/X11R6/include
      # Read from the CMakeSetup registry entries.  It is likely that
      # FLTK will have been recently built.
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild1]
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild2]
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild3]
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild4]
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild5]
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild6]
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild7]
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild8]
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild9]
      [HKEY_CURRENT_USER\\Software\\Kitware\\CMakeSetup\\Settings\\StartPath;WhereBuild10]
      # Look in standard Windows install locations.
      "C:/Program Files/FLTK/CMake"
    # Help the user find it if we cannot.
    DOC "The ${_FLTK_DIR_STRING}."
  )

  # If that wasn't successful then we'll look for other FLTK files
  FIND_PATH(FLTK_DIR
    NAMES
      include/FL/Fl.h
      include/FL/Fl.H
      bin/fltk-config
      bin/fluid
      lib/libfltk.a
      lib/libfltk.so
    PATHS
      # Consider environment variable first.
      $ENV{FLTK_DIR}
      # Look in places relative to the system executable search path.
      ${_FLTK_DIR_SEARCH}
      # Look in standard UNIX install locations.
      /usr
      /usr/local/lib/fltk
      /usr/lib/fltk
      /usr/local/include
      /usr/include
      /usr/local/fltk
      /usr/X11R6/include
      # Look in standard Windows install locations.
      "C:/Program Files/FLTK"
    # Help the user find it if we cannot.
    DOC "The ${_FLTK_DIR_STRING}."
  )

  UNSET(_FLTK_DIR_SEARCH_PATH)
  UNSET(_FLTK_CONFIG_SEARCH)
  UNSET(_FLTK_DIR_SEARCH)

ENDIF(NOT FLTK_DIR)

# Look for FLTK CMake configuration file or installed files and directories.
IF(FLTK_DIR)

  IF(EXISTS ${FLTK_DIR}/FLTKConfig.cmake OR EXISTS ${FLTK_DIR}/CMake/FLTKConfig.cmake)
    SET(FLTK_BUILT_WITH_CMAKE 1)
  ENDIF(EXISTS ${FLTK_DIR}/FLTKConfig.cmake OR EXISTS ${FLTK_DIR}/CMake/FLTKConfig.cmake)

  IF(FLTK_BUILT_WITH_CMAKE)

    SET(FLTK_FOUND 1)
    IF(EXISTS ${FLTK_DIR}/CMake/FLTKConfig.cmake)
      INCLUDE(${FLTK_DIR}/CMake/FLTKConfig.cmake)
    ELSE(EXISTS ${FLTK_DIR}/CMake/FLTKConfig.cmake)
      INCLUDE(${FLTK_DIR}/FLTKConfig.cmake)
    ENDIF(EXISTS ${FLTK_DIR}/CMake/FLTKConfig.cmake)

    IF(NOT FLTK_INCLUDE_DIR AND FLTK_INCLUDE_DIRS)
      SET(FLTK_INCLUDE_DIR "${FLTK_INCLUDE_DIRS}")
    ENDIF(NOT FLTK_INCLUDE_DIR AND FLTK_INCLUDE_DIRS)

    SET(FLTK_BASE_LIBRARY   fltk)
    SET(FLTK_GL_LIBRARY     fltk_gl)
    SET(FLTK_FORMS_LIBRARY  fltk_forms)
    SET(FLTK_IMAGES_LIBRARY fltk_images)

  ELSE(FLTK_BUILT_WITH_CMAKE)

    FIND_PROGRAM(FLTK_CONFIG_SCRIPT    NAMES fltk-config PATHS ${FLTK_DIR}/bin ${FLTK_DIR}/bin/Release ${FLTK_DIR}/bin/Debug ${FLTK_DIR}       NO_DEFAULT_PATH)
    FIND_PROGRAM(FLTK_FLUID_EXECUTABLE NAMES fluid       PATHS ${FLTK_DIR}/bin ${FLTK_DIR}/bin/Release ${FLTK_DIR}/bin/Debug ${FLTK_DIR}/fluid NO_DEFAULT_PATH)

    FIND_PATH(FLTK_INCLUDE_DIR NAMES FL/Fl.h FL/Fl.H PATHS ${FLTK_DIR}/include ${FLTK_DIR} NO_DEFAULT_PATH)

    FIND_LIBRARY(FLTK_BASE_LIBRARY   NAMES fltk       fltkd                   PATHS ${FLTK_DIR}/lib ${FLTK_DIR}/lib/Release ${FLTK_DIR}/lib/Debug NO_DEFAULT_PATH)
    FIND_LIBRARY(FLTK_GL_LIBRARY     NAMES fltkgl     fltkgld     fltk_gl     PATHS ${FLTK_DIR}/lib ${FLTK_DIR}/lib/Release ${FLTK_DIR}/lib/Debug NO_DEFAULT_PATH)
    FIND_LIBRARY(FLTK_FORMS_LIBRARY  NAMES fltkforms  fltkformsd  fltk_forms  PATHS ${FLTK_DIR}/lib ${FLTK_DIR}/lib/Release ${FLTK_DIR}/lib/Debug NO_DEFAULT_PATH)
    FIND_LIBRARY(FLTK_IMAGES_LIBRARY NAMES fltkimages fltkimagesd fltk_images PATHS ${FLTK_DIR}/lib ${FLTK_DIR}/lib/Release ${FLTK_DIR}/lib/Debug NO_DEFAULT_PATH)

    # extra libraries needed by the fltk_images library
    IF(UNIX AND FLTK_CONFIG_SCRIPT)
      EXEC_PROGRAM(${FLTK_CONFIG_SCRIPT} ARGS --use-images --ldflags OUTPUT_VARIABLE FLTK_IMAGES_LDFLAGS)
      SET(FLTK_LIBS_EXTRACT_REGEX ".*-lfltk_images (.*) -lfltk.*")
      IF(FLTK_IMAGES_LDFLAGS MATCHES "${FLTK_LIBS_EXTRACT_REGEX}")
        STRING(REGEX REPLACE "${FLTK_LIBS_EXTRACT_REGEX}" "\\1" FLTK_IMAGES_LIBRARIES "${FLTK_IMAGES_LDFLAGS}")
        STRING(REGEX REPLACE " +" ";" FLTK_IMAGES_LIBRARIES "${FLTK_IMAGES_LIBRARIES}")
        # The EXEC_PROGRAM will not be inherited into subdirectories from
        # the file that originally included this module.  Save the answer.
        SET(FLTK_IMAGES_LIBRARIES "${FLTK_IMAGES_LIBRARIES}" CACHE INTERNAL "Extra libraries for fltk_images library.")
      ENDIF(FLTK_IMAGES_LDFLAGS MATCHES "${FLTK_LIBS_EXTRACT_REGEX}")
    ENDIF(UNIX AND FLTK_CONFIG_SCRIPT)

  ENDIF(FLTK_BUILT_WITH_CMAKE)

ENDIF(FLTK_DIR)

# Determine whether all required files were found.
SET(FLTK_FOUND 1)
SET(_FLTK_MISSING_VARS)
FOREACH(_FLTK_VAR IN ITEMS FLTK_INCLUDE_DIR FLTK_BASE_LIBRARY FLTK_GL_LIBRARY FLTK_FORMS_LIBRARY FLTK_IMAGES_LIBRARY)
  IF(NOT ${_FLTK_VAR})
    SET(FLTK_FOUND 0)
    LIST(APPEND _FLTK_MISSING_VARS ${_FLTK_VAR})
  ENDIF(NOT ${_FLTK_VAR})
ENDFOREACH(_FLTK_VAR)

# Set uncached variables, in particular, FLTK_LIBRARIES list.
IF(FLTK_FOUND)

  SET(FLTK_LIBRARIES ${FLTK_IMAGES_LIBRARY} ${FLTK_IMAGES_LIBRARIES} ${FLTK_BASE_LIBRARY} ${FLTK_GL_LIBRARY} ${FLTK_FORMS_LIBRARY})
  IF(APPLE)
    SET(FLTK_LIBRARIES ${_FLTK_PLATFORM_DEPENDENT_LIBRARIES} ${FLTK_LIBRARIES})
  ELSE(APPLE)
    SET(FLTK_LIBRARIES ${FLTK_LIBRARIES} ${_FLTK_PLATFORM_DEPENDENT_LIBRARIES})
  ENDIF(APPLE)

  IF(UNIX AND NOT FLTK_BUILT_WITH_CMAKE)
    FIND_PACKAGE(X11 QUIET)
    IF (X11_Xft_LIB)
      SET(FLTK_LIBRARIES ${FLTK_LIBRARIES} ${X11_Xft_LIB})
    ENDIF (X11_Xft_LIB)
  ENDIF(UNIX AND NOT FLTK_BUILT_WITH_CMAKE)

  IF(FLTK_FLUID_EXECUTABLE)
    SET(FLTK_WRAP_UI 1)
  ELSE(FLTK_FLUID_EXECUTABLE)
    SET(FLTK_WRAP_UI 0)
  ENDIF(FLTK_FLUID_EXECUTABLE)

ENDIF(FLTK_FOUND)

# Mark cache entries set by this module as advanced (except of the FLTK_DIR).
FOREACH(_FLTK_VAR IN ITEMS
    FLTK_CONFIG_SCRIPT
    FLTK_FLUID_EXECUTABLE
    FLTK_INCLUDE_DIR
    FLTK_BASE_LIBRARY 
    FLTK_GL_LIBRARY
    FLTK_FORMS_LIBRARY
    FLTK_IMAGES_LIBRARY)
  # hide all advanced variables if set or FLTK_DIR not specified yet
  IF(${_FLTK_VAR} OR NOT FLTK_DIR)
    MARK_AS_ADVANCED(FORCE ${_FLTK_VAR})
  # otherwise, show advanced variable to the user so they are aware it has to be set
  ELSE(${_FLTK_VAR} OR NOT FLTK_DIR)
    MARK_AS_ADVANCED(CLEAR ${_FLTK_VAR})
  ENDIF(${_FLTK_VAR} OR NOT FLTK_DIR)
ENDFOREACH(_FLTK_VAR)
MARK_AS_ADVANCED(FLTK_LIBRARIES)

# Handle standard arguments.
IF(FLTK_FOUND)
  IF(NOT FLTK_FIND_QUIETLY)
    MESSAGE(STATUS "Looking for FLTK -- found")
  ENDIF(NOT FLTK_FIND_QUIETLY)
ELSE(FLTK_FOUND)
  IF(FLTK_FIND_REQUIRED)
    IF(FLTK_DIR)
      MESSAGE(FATAL_ERROR "Could NOT find FLTK! Please choose a proper FLTK_DIR by setting"
                          " it to the ${_FLTK_DIR_STRING} or set the following"
                          " variables:\n\t${_FLTK_MISSING_VARS}\n")
    ELSE(FLTK_DIR)
      MESSAGE(FATAL_ERROR "Could NOT find FLTK!\nPlease set FLTK_DIR to the ${_FLTK_DIR_STRING}.")
    ENDIF(FLTK_DIR)
  ELSE(FLTK_FIND_REQUIRED)
    IF(NOT FLTK_FIND_QUIETLY)
      MESSAGE( STATUS "Looking for FLTK -- not found" )
    ENDIF(NOT FLTK_FIND_QUIETLY)
  ENDIF(FLTK_FIND_REQUIRED)
ENDIF(FLTK_FOUND)

# Unset private variables.
UNSET(_FLTK_VAR)
UNSET(_FLTK_DIR_STRING)
UNSET(_FLTK_PLATFORM_DEPENDENT_LIBRARIES)
UNSET(_FLTK_MISSING_VARS)