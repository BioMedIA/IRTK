SET(COMMON_INCLUDES
../include/irtkCifstream.h
../include/irtkException.h
../include/irtkCofstream.h
../include/irtkObject.h
../include/irtkCommon.h
../include/irtkParallel.h
)

SET(COMMON_SRCS 
basename.cc 
dirname.cc 
irtkCommon.cc
irtkCifstream.cc
irtkException.cc 
irtkCofstream.cc 
irtkObject.cc
irtkParallel.cc
read.cc 
swap.cc
 weightedmedian.cc uncompress.cc)

ADD_LIBRARY(common++ ${COMMON_INCLUDES} ${COMMON_SRCS})

INSTALL_FILES(/include FILES ${COMMON_INCLUDES})

