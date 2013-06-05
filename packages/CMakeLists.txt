SUBDIRS(transformation segmentation registration registration2 rview applications)

IF (BUILD_CARDIAC)
  SUBDIRS(cardiac)
ENDIF (BUILD_CARDIAC)


