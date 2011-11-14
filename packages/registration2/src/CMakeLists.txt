SET(REGISTRATION2_INCLUDES
../include/irtkImageAffineRegistration2.h
../include/irtkImageFreeFormRegistration2.h
../include/irtkImageRegistration2.h
../include/irtkImageRigidRegistration2.h
../include/irtkDiscontinueFreeFormRegistration.h
)

SET(REGISTRATION2_SRCS 
irtkImageAffineRegistration2.cc
irtkImageFreeFormRegistration2.cc
irtkImageRegistration2.cc
irtkImageRigidRegistration2.cc
irtkDiscontinueFreeFormRegistration.cc
)

ADD_LIBRARY(registration2++ ${REGISTRATION2_INCLUDES} ${REGISTRATION2_SRCS})

INSTALL_FILES(/include FILES ${REGISTRATION2_INCLUDES})
