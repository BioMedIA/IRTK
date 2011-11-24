SET(CONTRIB_SRCS 
irtkArith.cc
irtkCityBlockDistanceTransform.cc
irtkEigenAnalysis.cc
irtkEuclideanDistanceTransform.cc
irtkMatrix3x3.cc
irtkVector3.cc)
ADD_LIBRARY(contrib++ ${CONTRIB_SRCS})
