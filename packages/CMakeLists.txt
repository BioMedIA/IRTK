SUBDIRS(transformation registration registration2 rview segmentation applications)

IF (BUILD_CARDIAC)
  SUBDIRS(cardiac)
ENDIF (BUILD_CARDIAC)

