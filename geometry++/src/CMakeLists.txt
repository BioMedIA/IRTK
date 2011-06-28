SET(GEOMETRY_INCLUDES
../include/irtkComplexFunction.h
../include/irtkGeometry.h
../include/irtkMatrix.h
../include/irtkNeighbourhoodOffsets.h
../include/irtkPoint.h
../include/irtkPointSet.h
../include/irtkQuaternion.h
../include/irtkScalarFunction.h
../include/irtkScalarGaussian.h
../include/irtkScalarGaussianDx.h
../include/irtkScalarGaussianDy.h
../include/irtkScalarGaussianDz.h
../include/irtkScalarGaussianDxDx.h
../include/irtkScalarGaussianDxDy.h
../include/irtkScalarGaussianDxDz.h
../include/irtkScalarGaussianDyDy.h
../include/irtkScalarGaussianDyDz.h
../include/irtkScalarGaussianDzDz.h
../include/irtkVector.h
../include/irtkVector3D.h
../include/irtkVTKFunctions.h)

SET(GEOMETRY_SRCS 
irtkComplexFunction.cc
irtkMatrix.cc
irtkNeighbourhoodOffsets.cc
irtkPoint.cc
irtkPointSet.cc
irtkQuaternion.cc
irtkScalarFunction.cc
irtkScalarGaussian.cc
irtkScalarGaussianDx.cc
irtkScalarGaussianDy.cc
irtkScalarGaussianDz.cc
irtkScalarGaussianDxDx.cc
irtkScalarGaussianDxDy.cc
irtkScalarGaussianDxDz.cc
irtkScalarGaussianDyDy.cc
irtkScalarGaussianDyDz.cc
irtkScalarGaussianDzDz.cc
irtkVector.cc
irtkVector3D.cc
irtkVTKFunctions.cc)

ADD_LIBRARY(geometry++ ${GEOMETRY_SRCS} ${GEOMETRY_INCLUDES})

INSTALL_FILES(/include FILES ${GEOMETRY_INCLUDES})
