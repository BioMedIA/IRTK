SUBDIRS(transformation registration segmentation rview applications)

IF (BUILD_CARDIAC)
  SUBDIRS(cardiac)
ENDIF (BUILD_CARDIAC)

