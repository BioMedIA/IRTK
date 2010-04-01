
IF (BUILD_CONDOR_EXE)
  #  SET (CMAKE_CXX_COMPILER "condor_compile ${CMAKE_CXX_COMPILER}")
  # This is a hack since the line above doesn't work: this local
  # scripts executes condor_compile with the option "c++". That means,
  # the variable CMAKE_CXX_COMPILER defined in ccmake is not taken
  # into account!
  SET (CMAKE_CXX_COMPILER "${IRTK_SOURCE_DIR}/condor_compile.sh")

  ADD_EXECUTABLE(ems.condor ems.cc)
  ADD_EXECUTABLE(emsbc.condor emsbc.cc)
  ADD_EXECUTABLE(bc2.condor bc2.cc)
  ADD_EXECUTABLE(ems_hard_segmentation.condor hard_segmentation.cc)


ELSE(BUILD_CONDOR_EXE)

  IF (BUILD_MPI_EXE)

  ELSE (BUILD_MPI_EXE)

    ADD_EXECUTABLE(bc bc.cc)
    ADD_EXECUTABLE(bc2 bc2.cc)
    ADD_EXECUTABLE(biasfield biasfield.cc)
    ADD_EXECUTABLE(biascorrect biascorrect.cc)
    ADD_EXECUTABLE(ems ems.cc)
    ADD_EXECUTABLE(emsbc emsbc.cc)
    ADD_EXECUTABLE(gmm gmm.cc)
    ADD_EXECUTABLE(gmm-par gmm-par.cc)
    ADD_EXECUTABLE(gmm-pv-par gmm-pv-par.cc)
    ADD_EXECUTABLE(ems_hard_segmentation hard_segmentation.cc)
    ADD_EXECUTABLE(enlarge_image enlarge_image.cc)

  ENDIF (BUILD_MPI_EXE)

INSTALL_TARGETS(/bin biasfield biascorrect ems emsbc)

ENDIF (BUILD_CONDOR_EXE)

