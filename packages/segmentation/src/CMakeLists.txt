SET(SEGMENTATION_INCLUDES
../include/irtkBiasCorrection.h
../include/irtkBiasField.h
../include/irtkBSplineBiasField.h
../include/irtkEMClassification.h
../include/irtkEMClassificationBiasCorrection.h
../include/irtkEMClassificationMultiComp.h
../include/irtkEMClassificationTemplateBiasCorrection.h
../include/irtkGaussian.h
../include/irtkIntensityMatching.h
../include/irtkMultiChannelImage.h
../include/irtkMultiImageGraphCut.h
../include/irtkMeanShift.h
../include/irtkProbabilisticAtlas.h
../include/irtkRician.h
../include/block.h
../include/graph.h
../include/instances.inc
../include/irtkImageGraphCut.h
../include/irtkSegmentationFunction.h)

SET(SEGMENTATION_SRCS
irtkBiasCorrection.cc
irtkBiasField.cc
irtkBSplineBiasField.cc
irtkEMClassification.cc
irtkEMClassificationBiasCorrection.cc
irtkEMClassificationMultiComp.cc
irtkEMClassificationTemplateBiasCorrection.cc
irtkGaussian.cc
irtkIntensityMatching.cc
irtkMultiChannelImage.cc
irtkMeanShift.cc
irtkProbabilisticAtlas.cc
irtkRician.cc
graph.cpp
irtkImageGraphCut.cc
irtkMultiImageGraphCut.cc
irtkSegmentationFunction.cc)

ADD_LIBRARY(segmentation++ ${SEGMENTATION_INCLUDES} ${SEGMENTATION_SRCS})
INSTALL_FILES(/include FILES ${SEGMENTATION_INCLUDES})
