#
# Build the documentation
#

INCLUDE (${CMAKE_ROOT}/Modules/Documentation.cmake OPTIONAL)

IF (BUILD_DOCUMENTATION)

  #
  # Configure the script and the doxyfile, then add target
  #
  CONFIGURE_FILE(
    ${IRTK_SOURCE_DIR}/documentation/doxyfile.in
    ${IRTK_BINARY_DIR}/documentation/doxyfile)

  ADD_CUSTOM_TARGET(doc ${DOXYGEN} doxyfile)

ENDIF (BUILD_DOCUMENTATION)

INSTALL_FILES(${IRTK_BINARY_DIR}/doc "\\.css$")
INSTALL_FILES(${IRTK_BINARY_DIR}/doc "\\.html$")
