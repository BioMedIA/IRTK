LINK_LIBRARIES(cardiac++ registration++ segmentation++ transformation++ contrib++ image++ geometry++ common++ contrib++ image++ geometry++ common++) 
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

ELSE(BUILD_CONDOR_EXE)

  IF (BUILD_MPI_EXE)

  ELSE (BUILD_MPI_EXE)
    ADD_EXECUTABLE(bep bep.cc)
    ADD_EXECUTABLE(cardiacdof2vtk cardiacdof2vtk.cc)
    ADD_EXECUTABLE(cardiacdotspacing cardiacdotspacing.cc)
    ADD_EXECUTABLE(cardiacsurfacevolume cardiacsurfacevolume.cc)
    ADD_EXECUTABLE(comparebep comparebep.cc)
    ADD_EXECUTABLE(ejecfrac ejecfrac.cc)
    ADD_EXECUTABLE(late late.cc)
    ADD_EXECUTABLE(motiontrackcardiac motiontrackcardiac.cc)
    ADD_EXECUTABLE(spatial_correct spatial_correct.cc)
    ADD_EXECUTABLE(tracktag tracktag.cc)

  ENDIF (BUILD_MPI_EXE)

INSTALL_TARGETS(/bin bep cardiacdof2vtk cardiacdotspacing cardiacsurfacevolume comparebep ejecfrac late motiontrackcardiac spatial_correct tracktag)

ENDIF (BUILD_CONDOR_EXE)

