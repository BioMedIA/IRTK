SET(REGISTRATION2_INCLUDES
../include/irtkImageAffineRegistration2.h
../include/irtkImageFreeFormRegistration2.h
../include/irtkImageRegistration2.h
../include/irtkImageRigidRegistration2.h
../include/irtkMultipleImageRegistration2.h
../include/irtkMultipleImageRigidRegistration2.h
../include/irtkMultipleImageFreeFormRegistration2.h
../include/irtkSparseFreeFormRegistration.h
../include/irtkTemporalImageRegistration.h
../include/irtkImageTFFDRegistration.h
)

SET(REGISTRATION2_SRCS 
irtkImageAffineRegistration2.cc
irtkImageFreeFormRegistration2.cc
irtkMultipleImageRegistration2.cc
irtkMultipleImageRigidRegistration2.cc
irtkMultipleImageFreeFormRegistration2
irtkImageRegistration2.cc
irtkImageRigidRegistration2.cc
irtkSparseFreeFormRegistration.cc
irtkTemporalImageRegistration.cc
irtkImageTFFDRegistration.cc
)

ADD_LIBRARY(registration2++ ${REGISTRATION2_INCLUDES} ${REGISTRATION2_SRCS})

INSTALL_FILES(/include FILES ${REGISTRATION2_INCLUDES})
