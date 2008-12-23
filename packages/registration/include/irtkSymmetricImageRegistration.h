/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSYMMETRICIMAGEREGISTRATION_H

#define _IRTKSYMMETRICIMAGEREGISTRATION_H

#include <irtkSymmetricOptimizer.h>

/**
 * Generic for symmetric image registration based on voxel similarity measures.
 *
 * This class implements a registration filter which takes two input images
 * and calculates two transformation. The first maps the second image (denoted
 * as source image) into the coordinate system of the first image (denoted as
 * target image).  The second transformation is defined in the opposite direction.
 * This is the abstract base class which defines a common interface for arbitrary
 * registration filters. Each derived class has to implement all abstract member
 * functions.
 *
 */

class irtkSymmetricImageRegistration : public irtkRegistration
{

  /// Interface to input file stream
  friend istream& operator>> (istream&, irtkSymmetricImageRegistration*);

  /// Interface to output file stream
  friend ostream& operator<< (ostream&, const irtkSymmetricImageRegistration*);

protected:

  /** First input image. This image is denoted as target image and its
   *  coordinate system defines the frame of reference for the registration.
   */
  irtkGreyImage *_target;

  /** Second input image. This image is denoted as source image. The goal of
   *  the registration is to find the transformation which maps the source
   *  image into the coordinate system of the target image.
   */
  irtkGreyImage *_source;

  /// Outputs
  irtkTransformation *_transformation1;
  irtkTransformation *_transformation2;

  /// Histogram
  irtkSimilarityMetric *_metric1;
  irtkSimilarityMetric *_metric2;

  /// Interpolator
  irtkInterpolateImageFunction<irtkGreyPixel> *_interpolator1;
  irtkInterpolateImageFunction<irtkGreyPixel> *_interpolator2;

  /// Optimizer
  irtkSymmetricOptimizer *_optimizer;

  /// Blurring of target image (in mm)
  double _TargetBlurring[MAX_NO_RESOLUTIONS];

  /// Resolution of target image (in mm)
  double _TargetResolution[MAX_NO_RESOLUTIONS][3];

  /// Blurring of source image (in mm)
  double _SourceBlurring[MAX_NO_RESOLUTIONS];

  /// Resolution of source image (in mm)
  double _SourceResolution[MAX_NO_RESOLUTIONS][3];

  /// Number of step sizes
  int    _NumberOfSteps[MAX_NO_RESOLUTIONS];

  /// Length of steps
  double _LengthOfSteps[MAX_NO_RESOLUTIONS];

  /// Max. number of iterations per step size
  int    _NumberOfIterations[MAX_NO_RESOLUTIONS];

  /// Padding value of target image
  short  _Padding;

  /// Number of levels of multiresolution pyramid
  int    _NumberOfLevels;

  /// Max. number of bins for histogram
  int    _NumberOfBins;

  /// Similarity measure for registration
  irtkSimilarityMeasure  _SimilarityMeasure;

  /// Optimization method for registration
  irtkOptimizationMethod _OptimizationMethod;

  /// Interpolation mode to use during resampling and registration
  irtkInterpolationMode _InterpolationMode;

  /// Convergence parameter for optimization based on change in similarity.
  double _Epsilon;

  /// Convergence parameter for optimization based on change in the transformation.
  double _Delta;

  /// Debugging flag
  int    _DebugFlag;

  /// Source image domain which can be interpolated fast
  double _source_x1, _source_y1, _source_z1;
  double _source_x2, _source_y2, _source_z2;

  /// Target image domain which can be interpolated fast
  double _target_x1, _target_y1, _target_z1;
  double _target_x2, _target_y2, _target_z2;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

  /// Initial set up for the registration at a multiresolution level
  virtual void Initialize(int);

  /// Final set up for the registration at a multiresolution level
  virtual void Finalize(int);

public:

  /// Classification
  irtkEMClassification *classification;

  /// Constructor
  irtkSymmetricImageRegistration();

  /// Destructor
  virtual ~irtkSymmetricImageRegistration();

  /// Sets input for the registration filter
  virtual void SetInput (irtkGreyImage *, irtkGreyImage *);

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *, irtkTransformation *) = 0;

  /// Runs the registration filter
  virtual void Run();

  /** Evaluates the similarity metric. This function evaluates the similarity
   *  metric of the registration by looping over the target image and
   *  interpolating the transformed source image while filling the joint
   *  histogram. This function returns the value of the similarity measure
   *  using Similarity().
   */
  virtual double Evaluate() = 0;

  /** Evaluates the gradient of the similarity metric. This function
   *  evaluates the gradient of the similarity metric of the registration
   *  by looping over the target image and interpolating the transformed
   *  source image while filling the joint histogram. The partial derivatives
   *  are approximated using a finite difference scheme. The step size for the
   *  finite difference scheme is passed as a parameter to the
   *  function. The function returns the norm of the gradient vector as
   *  well as the gradient vector containing the partial derivatives.
   */
  virtual double EvaluateGradient(float, float *);

  /// Returns the name of the class
  virtual const char *NameOfClass() = 0;

  /// Prints debugging messages if debugging is enabled
  virtual void Debug(string);

  /// Prints information about the progress of the registration
  virtual void Print() = 0;

  /// Guess parameters
  virtual void GuessParameter() = 0;

  /// Read registration parameters from file
  virtual void Read (char *);

  /// Parse parameter line
  virtual Bool Read(char *, char *, int &);

  /// Write registration parameters to file
  virtual void Write(char *);

  /// Write parameters to stream
  virtual void Write(ostream &);

  // Access parameters
  virtual SetMacro(DebugFlag, int);
  virtual GetMacro(DebugFlag, int);
  virtual SetMacro(Padding, int);
  virtual GetMacro(Padding, int);
  virtual SetMacro(OptimizationMethod, irtkOptimizationMethod);
  virtual GetMacro(OptimizationMethod, irtkOptimizationMethod);

};

inline void irtkSymmetricImageRegistration::SetInput(irtkGreyImage *target, irtkGreyImage *source)
{
  _target = target;
  _source = source;
}

inline void irtkSymmetricImageRegistration::Debug(string message)
{
  if (_DebugFlag == True) cout << message << endl;
}

#include <irtkSymmetricImageFreeFormRegistration.h>

#endif
