PROJECT(RVIEW_SRC)
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
