PROJECT(RVIEW_SRC)
SET(RVIEW_INCLUDES
../include/irtkColor.h
../include/irtkColorRGBA.h
../include/irtkContour.h
../include/irtkLookupTable.h
../include/irtkRView.h
../include/irtkRViewConfig.h
../include/irtkViewer.h
../include/irtkHistogramWindow.h
../include/irtkSegment.h
../include/irtkSegmentTable.h
../include/irtkVoxelContour.h)

SET(RVIEW_SRCS
irtkColor.cc
irtkColorRGBA.cc
irtkLookupTable.cc
irtkRView.cc
irtkRViewConfig.cc
irtkViewer.cc
irtkHistogramWindow.cc
irtkSegment.cc
irtkSegmentTable.cc
irtkVoxelContour.cc)

ADD_LIBRARY(rview++ ${RVIEW_SRCS})

INSTALL_FILES(/include FILES ${RVIEW_INCLUDES})