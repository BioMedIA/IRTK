SET(TRANSFORMATION_INCLUDES
../include/irtkAffineTransformation.h
../include/irtkBSplineFreeFormTransformation3D.h
../include/irtkBSplineFreeFormTransformation4D.h
../include/irtkEigenFreeFormTransformation.h
../include/irtkFreeFormTransformation.h
../include/irtkFreeFormTransformation3D.h
../include/irtkFreeFormTransformation4D.h
../include/irtkFluidFreeFormTransformation.h
../include/irtkHomogeneousTransformation.h
../include/irtkHomogeneousTransformationIterator.h
../include/irtkImageHomogeneousTransformation.h
../include/irtkImageFastFourierTransform.h
../include/irtkImageTransformation.h
../include/irtkImageTransformation2.h
../include/irtkLinearFreeFormTransformation.h
../include/irtkMultiLevelFreeFormTransformation.h
../include/irtkQuaternionTransformation.h
../include/irtkRigidTransformation.h
../include/irtkTransformation.h
../include/irtkGradVecFlow.h
../include/newt2.h)

SET(TRANSFORMATION_SRCS
irtkAffineTransformation.cc
irtkHomogeneousTransformation.cc
irtkImageHomogeneousTransformation.cc
irtkImageTransformation.cc
irtkImageTransformation2.cc
irtkImageFastFourierTransform.cc
irtkEigenFreeFormTransformation.cc
irtkBSplineFreeFormTransformation3D.cc
irtkBSplineFreeFormTransformation4D.cc
irtkLinearFreeFormTransformation.cc
irtkFreeFormTransformation.cc
irtkFreeFormTransformation3D.cc
irtkFreeFormTransformation4D.cc
irtkFluidFreeFormTransformation.cc
irtkMultiLevelFreeFormTransformation.cc
irtkQuaternionTransformation.cc
irtkRigidTransformation.cc
irtkTransformation.cc
irtkGradVecFlow.cc
newt2.cc)

ADD_LIBRARY(transformation++ ${TRANSFORMATION_SRCS} ${TRANSFORMATION_INCLUDES})

INSTALL_FILES(/include FILES ${TRANSFORMATION_INCLUDES})
