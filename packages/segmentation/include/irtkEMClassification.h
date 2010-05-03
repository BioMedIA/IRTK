/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKEMCLASSIFICATION_H

#define _IRTKEMCLASSIFICATION_H

#include <irtkImage.h>

#include <irtkProbabilisticAtlas.h>
#include <irtkGaussian.h>

/*

Expectation Maximisation Algorithm

*/

class irtkEMClassification : public irtkObject
{

protected:

  /// Input image
  irtkRealImage _input;

  /// Input mask
  irtkRealImage _mask;

  /// Brain tissue
  irtkRealImage _brain;

  /// weights image
  irtkRealImage _weights;

  /// weights image
  irtkRealImage _weightsB;

  /// weights image
  irtkRealImage _weightsR;

  /// image estimate
  irtkRealImage _estimate;

  /// Posterior probability maps -  segmentation
  irtkProbabilisticAtlas _output;

  /// Partial volume soft segmentation
  irtkProbabilisticAtlas _pv_output;

  /// Probability maps (atlas)
  irtkProbabilisticAtlas _atlas;

  /// image segmentation
  irtkRealImage _segmentation;

  /// distance map
  irtkRealImage _distance;

  ///Subcortical mask
  irtkGreyImage _smask;




  /// Number of tissues
  int _number_of_tissues;

  /// Number of voxels
  int _number_of_voxels;

  /// Gaussian distribution parameters for each tissue type
  double *_mi;
  double *_sigma;
  /// mixing coefficients for GMM
  double *_c;

  irtkGaussian *_G;

  /// Likelihood
  double _f;

  /// Padding value
  irtkRealPixel _padding;

  int _background_tissue;

  bool _debug;


public:

  /// Estimates posterior probabilities
  virtual void EStep();

  /// Estimates parameters
  virtual void MStep();

  /// Estimates posterior probabilities
  virtual void EStepGMM(bool uniform_prior = false);

  /// Estimates parameters in GMM
  virtual void MStepGMM(bool uniform_prior = false);

  /// Estimates parameters in GMM with equal variance
  virtual void MStepVarGMM(bool uniform_prior = false);

  /// Remove partial volume misclassifications in neonatal segmentation
  virtual bool PVStep(int wm1Label, int wm2Label, int cortexLabel, int csfLabel, int backgroundLabel, double lcc_treshold = 0.5, bool first = false, irtkGreyImage * eyes = NULL);

  /// Computes log likelihood for current parameters
  virtual double LogLikelihood();

  /// Computes log likelihood for GMM for current parameters
  virtual double LogLikelihoodGMM();

public:
  /// Print parameters
  virtual void Print();

  /// Print parameters
  virtual void PrintGMM();

  /// Calculate weights
  virtual void WStep();

  /// Empty constructor
  irtkEMClassification();

  /// Constructor when adding background
  irtkEMClassification(int noTissues, irtkRealImage **atlas, irtkRealImage *background);

  /// Constructor without adding background
  irtkEMClassification(int noTissues, irtkRealImage **atlas);

  /// Destructor
  virtual ~irtkEMClassification();

  /// Initialize filter
  virtual void Initialise();

  /// Initialize filter
  virtual void InitialiseGMM();

  /// Initialize posteriors if GMM parameters were given first
  void InitialisePosteriors();

  /// Initialize atlas if no probability maps were given
  void InitialiseAtlas();

  /// Initialize PV segmentation
  void InitialisePVSegmentation();

  /// Initialize Gaussian Mixture Parameters and calculate initial posteriors
  void InitialiseGMMParameters(int n, double *m, double *s, double *c);

  /// Set image
  virtual void SetInput(const irtkRealImage &);

  /// Set padding value
  virtual void SetPadding(irtkRealPixel);

  /// Execute one iteration and return log likelihood
  virtual double Iterate(int iteration);

  /// Execute one iteration and return log likelihood
  virtual double IterateGMM(int iteration, bool equal_var = false, bool uniform_prior = false);

  /// Construct segmentation based on current posterior probabilities
  virtual void ConstructSegmentation(irtkRealImage &);
  /// Construct segmentation based on current posterior probabilities
  virtual void ConstructSegmentation();
  /// Construct segmentation based on current PV segmentation
  virtual void ConstructSegmentationFromPV();
  /// Construct segmentation based on current posterior probabilities
  virtual void ConstructSegmentationWithPadding(irtkRealImage &);
  /// Construct segmentation based on current posterior probabilities
  virtual void ConstructSegmentationNoBG(irtkRealImage &);
  /// Construct PV soft segmentation based on current segmentation and posterior probabilities
  virtual void ConstructPVSegmentation();
  /// Construct PV segmentation with specific application for neonatal brain
  void ConstructSegmentationBrainNonBrain(int wm1Label, int wm2Label, int cortexLabel, int csfLabel, int backgroundLabel, irtkGreyImage* smask);



  ///Write probability map into a file
  void WriteProbMap(int i, char *filename);
  ///Write Gaussian parameters into a file
  void WriteGaussianParameters(char *file_name);
  ///Write image estimate
  void WriteEstimate(char *filename);
  ///Write weights
  void WriteWeights(char *filename);
  ///Write input
  void WriteInput(char *filename);
  ///Write segmentation
  void WriteSegmentation(char *filename);
  ///Write hard PV segmentation
  void WritePVSegmentation();
  ///Write soft PV segmentation
  void WritePVProbMap(int i, char *filename);


  ///Returns log likelihood for given intensity value
  double PointLogLikelihoodGMM(double x);
  virtual double PointLogLikelihoodGMM(double x, double y);
  virtual double PointLogLikelihoodGMMnomatch(double x, double y);
  ///Initialize gaussians with current GMM parameters
  void GInit();
  void CreateMask();
  void SetMask(irtkRealImage &mask);
  virtual double MatchIntensity(double) {
    return 0;
  }
  virtual void BrainmaskInput();
  inline double LLGMM();

  void InitialiseGMMParameters3();
  void InitialiseGMMParameters(int n);
  void UniformPrior();
  void DistanceTransform( int label1, int label2);

  virtual void SetDebugFlag(bool debug);

};

inline void irtkEMClassification::SetDebugFlag(bool debug)
{
  _debug = debug;
}

inline void irtkEMClassification::SetPadding(irtkRealPixel padding)
{
  _padding = padding;
}

inline void irtkEMClassification::WriteEstimate(char *filename)
{
  _estimate.Write(filename);
}

inline void irtkEMClassification::WriteInput(char *filename)
{
  _input.Write(filename);
}

inline void irtkEMClassification::WriteSegmentation(char *filename)
{
  _segmentation.Write(filename);
}

inline double irtkEMClassification::LLGMM()
{
  return _f;
}


#endif
