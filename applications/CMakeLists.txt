LINK_LIBRARIES(contrib++ image++ geometry++ common++) 
IF (USE_VXL)
  LINK_LIBRARIES(irtkvnl irtkvcl irtkvnl_algo irtknetlib)
ELSE (USE_VXL)
  LINK_LIBRARIES(recipes)
ENDIF (USE_VXL)



IF (BUILD_CONDOR_EXE)
  #  SET (CMAKE_CXX_COMPILER "condor_compile ${CMAKE_CXX_COMPILER}")
  # This is a hack since the line above doesn't work: this local 
  # scripts executes condor_compile with the option "c++". That means,
  # the variable CMAKE_CXX_COMPILER defined in ccmake is not taken 
  # into account!
  SET (CMAKE_CXX_COMPILER "${IRTK_SOURCE_DIR}/condor_compile.sh")

  ADD_EXECUTABLE(convert.condor convert.cc)
  ADD_EXECUTABLE(dilation.condor dilation.cc)
  ADD_EXECUTABLE(erosion.condor erosion.cc)
  ADD_EXECUTABLE(info.condor info.cc)
  ADD_EXECUTABLE(noise.condor noise.cc)
  ADD_EXECUTABLE(padding.condor padding.cc)
  ADD_EXECUTABLE(region.condor region.cc)

ELSE(BUILD_CONDOR_EXE)

  IF (BUILD_MPI_EXE)

  ELSE (BUILD_MPI_EXE)

    ADD_EXECUTABLE(combineLabels combineLabels.cc)
    ADD_EXECUTABLE(convert convert.cc)
    ADD_EXECUTABLE(dilation dilation.cc)
    ADD_EXECUTABLE(erosion erosion.cc)
    ADD_EXECUTABLE(info info.cc)
    ADD_EXECUTABLE(noise noise.cc)
    ADD_EXECUTABLE(padding padding.cc)
    ADD_EXECUTABLE(region region.cc)

    ADD_EXECUTABLE(binarize binarize.cc)
    ADD_EXECUTABLE(blur blur.cc)
    ADD_EXECUTABLE(closing closing.cc)
    ADD_EXECUTABLE(dicom2image dicom2image.cc)
    ADD_EXECUTABLE(dmap dmap.cc)
    ADD_EXECUTABLE(edgedetect edgedetect.cc)
    ADD_EXECUTABLE(flip flip.cc)
    ADD_EXECUTABLE(giplmagic giplmagic.cc)
    ADD_EXECUTABLE(headertool headertool.cc)
    ADD_EXECUTABLE(labelStats labelStats.cc)
    ADD_EXECUTABLE(lcc lcc.cc)
    ADD_EXECUTABLE(lcc_iterative lcc_iterative.cc)
    ADD_EXECUTABLE(makevolume makevolume.cc)
    ADD_EXECUTABLE(makesequence makesequence.cc)
    ADD_EXECUTABLE(mcubes mcubes.cc)
    ADD_EXECUTABLE(msample msample.cc)
    ADD_EXECUTABLE(opening opening.cc)
    ADD_EXECUTABLE(overlap overlap.cc)
    ADD_EXECUTABLE(raw2image raw2image.cc)
    ADD_EXECUTABLE(reflect reflect.cc)
    ADD_EXECUTABLE(resample resample.cc)
    ADD_EXECUTABLE(rescale rescale.cc)
    ADD_EXECUTABLE(rescale_real rescale_real.cc)
    ADD_EXECUTABLE(threshold threshold.cc)
    ADD_EXECUTABLE(voxelsize voxelsize.cc)
    ADD_EXECUTABLE(histogram histogram.cc)
    ADD_EXECUTABLE(sbig  sbig.cc)	
    ADD_EXECUTABLE(anisoDiff AnisoDiff.cc)
    
  ENDIF (BUILD_MPI_EXE)

INSTALL_TARGETS(/bin binarize blur closing convert dilation
dicom2image dmap edgedetect erosion flip giplmagic headertool info lcc
lcc_iterative makesequence makevolume mcubes msample noise opening 
overlap padding raw2image reflect region resample threshold voxelsize 
histogram sbig anisoDiff)

ENDIF (BUILD_CONDOR_EXE)

