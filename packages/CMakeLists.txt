#SUBDIRS(transformation registration registration2 rview segmentation applications applications_stefan)
SUBDIRS(transformation registration registration2 rview applications applications_stefan)

IF (BUILD_CARDIAC)
  SUBDIRS(cardiac)
ENDIF (BUILD_CARDIAC)

