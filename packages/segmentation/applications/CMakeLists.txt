
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

    ADD_EXECUTABLE(atrophy_graphcut4D atrophy_graphcut4D.cc)
    ADD_EXECUTABLE(biasfield biasfield.cc)
    ADD_EXECUTABLE(biascorrect biascorrect.cc)
    ADD_EXECUTABLE(createsamples createsamples.cc)
    ADD_EXECUTABLE(detect detect.cc)
    ADD_EXECUTABLE(ems ems.cc)
    ADD_EXECUTABLE(emsmc emsmc.cc)
    ADD_EXECUTABLE(emsbc emsbc.cc)
    ADD_EXECUTABLE(gmm gmm.cc)
    ADD_EXECUTABLE(gmm-par gmm-par.cc)
    ADD_EXECUTABLE(gmm-pv-par gmm-pv-par.cc)
    ADD_EXECUTABLE(graphcut graphcut.cc)
    ADD_EXECUTABLE(ems_hard_segmentation hard_segmentation.cc)
    ADD_EXECUTABLE(enlarge_image enlarge_image.cc)
    ADD_EXECUTABLE(multiimage_graphcut multiimage_graphcut.cc)
    ADD_EXECUTABLE(patchBasedSegmentation patchBasedSegmentation.cc)
    ADD_EXECUTABLE(patchBasedSegmentation_native patchBasedSegmentation_native.cc)
    ADD_EXECUTABLE(patchmatch patchmatch.cc)
    ADD_EXECUTABLE(patchmatchSegmentation patchmatchSegmentation.cc)
    ADD_EXECUTABLE(seg_uncertainty seg_uncertainty.cc)
    ADD_EXECUTABLE(reconstruction reconstruction.cc)
    ADD_EXECUTABLE(emSeg2ndOrderMRF emSeg2ndOrderMRF.cc)

  ENDIF (BUILD_MPI_EXE)

INSTALL_TARGETS(/bin emSeg2ndOrderMRF atrophy_graphcut4D biasfield biascorrect createsamples detect ems emsbc emsmc graphcut multiimage_graphcut patchBasedSegmentation patchBasedSegmentation_native patchmatch patchmatchSegmentation seg_uncertainty)

ENDIF (BUILD_CONDOR_EXE)

