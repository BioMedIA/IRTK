SET(SEGMENTATION_INCLUDES
../include/irtkBiasCorrection.h
../include/irtkBiasField.h
../include/irtkBSplineBiasField.h
../include/irtkBSplineReconstruction.h
../include/irtkEMClassification.h
../include/irtkEMClassificationBiasCorrection.h
../include/irtkEMClassificationBiasCorrectionfMRI.h
../include/irtkEMClassificationMultiComp.h
../include/irtkGaussian.h
../include/irtkIntensityMatching.h
../include/irtkMultiChannelImage.h
../include/irtkMultiImageGraphCut.h
../include/irtkMeanShift.h
../include/irtkProbabilisticAtlas.h
../include/irtkPatchMatch.h
../include/irtkMAPatchMatch.h
../include/irtkMAPatchMatchSegmentation.h
../include/irtkMAPatchMatchSuperResolution.h
../include/irtkReconstruction.h
../include/irtkRician.h
../include/irtkSubcorticalSegmentation_4D.h
../../../external/gco-v3.0/block.h
../../../external/gco-v3.0/graph.h
../../../external/gco-v3.0/energy.h
../../../external/gco-v3.0/GCoptimization.h
../../../external/gco-v3.0/LinkedBlockList.h
../include/irtkCRF.h
../include/irtkGraphCutSegmentation_4D.h
../include/irtkImageGraphCut.h
../include/irtkSegmentationFunction.h
../include/irtkPatchBasedSegmentation.h
../include/irtkEMClassification2ndOrderMRF.h
../include/irtkBiasCorrectionMask.h
../include/irtkPolynomialBiasField.h)

SET(SEGMENTATION_SRCS
irtkBiasCorrection.cc
irtkBiasField.cc
irtkBSplineBiasField.cc
irtkBSplineReconstruction.cc
irtkEMClassification.cc
irtkEMClassificationBiasCorrection.cc
irtkEMClassificationBiasCorrectionfMRI.cc
irtkEMClassificationMultiComp.cc
irtkGaussian.cc
irtkIntensityMatching.cc
irtkMultiChannelImage.cc
irtkMeanShift.cc
irtkProbabilisticAtlas.cc
irtkPatchMatch.cc
irtkMAPatchMatch.cc
irtkMAPatchMatchSegmentation.cc
irtkMAPatchMatchSuperResolution.cc
irtkRician.cc
irtkReconstruction.cc
irtkSubcorticalSegmentation_4D.cc
irtkCRF.cc
../../../external/gco-v3.0/graph.cpp
../../../external/gco-v3.0/GCoptimization.cpp
../../../external/gco-v3.0/LinkedBlockList.cpp
../../../external/gco-v3.0/maxflow.cpp
irtkGraphCutSegmentation_4D.cc
irtkImageGraphCut.cc
irtkMultiImageGraphCut.cc
irtkPatchBasedSegmentation.cc
irtkSegmentationFunction.cc
irtkEMClassification2ndOrderMRF.cc
irtkBiasCorrectionMask.cc
irtkPolynomialBiasField.cc)


IF (BUILD_GPU_SUPPORT_WITH_CUDA)
#TODO 
#full integration into the IRTK build environment does not work yet 
#-> build separately for now
#add_subdirectory(irtkReconstructionCuda)
ENDIF(BUILD_GPU_SUPPORT_WITH_CUDA)

ADD_LIBRARY(segmentation++ ${SEGMENTATION_INCLUDES} ${SEGMENTATION_SRCS})
INSTALL_FILES(/include FILES ${SEGMENTATION_INCLUDES})

