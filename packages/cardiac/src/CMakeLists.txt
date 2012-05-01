SET(CARDIAC_INCLUDES
../include/irtkCardiac.h
../include/irtkBep.h
../include/irtkCardiac3DImageFreeFormRegistration.h
../include/irtkCardiacSpatialCorrection.h
../include/irtkTagFunction.h
)

SET(CARDIAC_SRCS
irtkBep.cc
irtkCardiac3DImageFreeFormRegistration.cc
irtkCardiacSpatialCorrection.cc
irtkTagFunction.cc
)

ADD_LIBRARY(cardiac++ ${CARDIAC_INCLUDES} ${CARDIAC_SRCS})
INSTALL_FILES(/include FILES ${CARDIAC_INCLUDES})
