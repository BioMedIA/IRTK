SET(GEOMETRY_INCLUDES
../include/irtkComplexFunction.h
../include/irtkGeometry.h
../include/irtkMatrix.h
../include/irtkPoint.h
../include/irtkPointSet.h
../include/irtkQuaternion.h
../include/irtkScalarFunction.h
../include/irtkScalarGaussian.h
../include/irtkVector.h
../include/irtkVector3D.h)

SET(GEOMETRY_SRCS 
irtkComplexFunction.cc
irtkMatrix.cc
irtkPoint.cc
irtkPointSet.cc
irtkQuaternion.cc
irtkScalarFunction.cc
irtkScalarGaussian.cc
irtkVector.cc
irtkVector3D.cc)

ADD_LIBRARY(geometry++ ${GEOMETRY_SRCS} ${GEOMETRY_INCLUDES})

INSTALL_FILES(/include FILES ${GEOMETRY_INCLUDES})
