/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkCardiacSpatialCorrection.h 293 2011-02-03 21:07:57Z ws207 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2011-02-03 21:07:57 +0000 (Thu, 03 Feb 2011) $
  Version   : $Revision: 293 $
  Changes   : $Author: ws207 $

=========================================================================*/

#ifndef _IRTKCARDIACSPATIALCORRECTION_H

#define _IRTKCARDIACSPATIALCORRECTION_H

/**
 * Generic for multiple transformation registration based on voxel similarity measures.
 *
 * This class implements a registration filter which takes two sets of input
 * images and calculates the transformation which maps the second set of images
 * (denote as source image) into the coordinate system of the first image set
 * (denotes astarget image).  This is the abstract base class which defines a
 * common interface for arbitrary registration filters. Each derived class has to
 * implement all abstract member functions.
 *
 */

class irtkCardiacSpatialCorrection : public irtkRegistration
{

  /// Interface to input file stream
  friend istream& operator>> (istream&, irtkCardiacSpatialCorrection*);

  /// Interface to output file stream
  friend ostream& operator<< (ostream&, const irtkCardiacSpatialCorrection*);

protected:

  /** First set of input image. This image is denoted as target image and its
   *  coordinate system defines the frame of reference for the registration.
   */
  irtkGreyImage *_target;

  /** Second input image. This image is denoted as source image. The goal of
   *  the registration is to find the transformation which maps the source
   *  image into the coordinate system of the target image.
   */
  irtkGreyImage *_source;

  /// Landmark regulation Input
  irtkPointSet *_ptarget;

  /// Landmark regulation Input
  irtkPointSet *_psource;

  /// Output
  irtkTransformation **_transformation;

  /// Initial transformation
  irtkTransformation *_origintransformation;

  /// Histogram
  irtkSimilarityMetric *_metric;

  /// Internal Histogram
  irtkSimilarityMetric *_intmetric;

  /// Interpolator
  irtkInterpolateImageFunction *_interpolator;

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
  short  _TargetPadding;

  /// Number of levels of multiresolution pyramid
  int    _NumberOfLevels;

  /// Max. number of bins for histogram
  int    _NumberOfBins;

  /// Similarity measure for registration
  irtkSimilarityMeasure  _SimilarityMeasure;

  /// Interpolation mode to use during resampling and registration
  irtkInterpolationMode _InterpolationMode;

  /// Convergence parameter for optimization based on change in similarity.
  double _Epsilon;

  /// Inter slice regulation term
  double _IntLambda;

  /// transformations regulation term
  double _ConstrainLambda;

  /// Convergence parameter for optimization based on change in the transformation.
  double _Delta[MAX_NO_RESOLUTIONS];

  /// Debugging flag
  int    _DebugFlag;

  /// Source image domain which can be interpolated fast
  double _source_x1, _source_y1, _source_z1;
  double _source_x2, _source_y2, _source_z2;

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
  irtkCardiacSpatialCorrection();

  /// Destructor
  virtual ~irtkCardiacSpatialCorrection();

  /// Sets input for the registration filter
  virtual void SetInput (irtkGreyImage *, irtkGreyImage *);

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation **);

  /// Runs the registration filter
  virtual void Run();

  /** Evaluates the similarity metric. This function evaluates the similarity
   *  metric of the registration by looping over the target image and
   *  interpolating the transformed source image while filling the joint
   *  histogram. This function returns the value of the similarity measure
   *  using Similarity().
   */
  virtual double Evaluate();

  /** Evaluates the similarity metric with in the image.
   */
  virtual double InternalSimilarity();

  /** Constrain the transformation degree.
   */
  virtual double ConstrainTransformation();

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
  virtual const char *NameOfClass();

  /// Prints debugging messages if debugging is enabled
  virtual void Debug(string);

  /// Read registration parameters from file
  virtual void Read (char *);

  /// Parse parameter line
  virtual bool Read(char *, char *, int &);

  /// Write registration parameters to file
  virtual void Write(char *);

  /// Write parameters to stream
  virtual void Write(ostream &);

  // Access parameters
  virtual SetMacro(DebugFlag, int);
  virtual GetMacro(DebugFlag, int);
  virtual SetMacro(TargetPadding, int);
  virtual GetMacro(TargetPadding, int);

  /// Guess parameter
  virtual void GuessParameter();

  /// Run Optimization
  virtual void OptimizerRun(double&,double&,double&);

  /// Run Optimization step
  virtual double RunStep(double&);

};

inline const char *irtkCardiacSpatialCorrection::NameOfClass()
{
	return "irtkCardiacSpatialCorrection";
}

inline void irtkCardiacSpatialCorrection::SetInput(irtkGreyImage *target, irtkGreyImage *source)
{
  _target = target;
  _source = source;
}

inline void irtkCardiacSpatialCorrection::SetOutput(irtkTransformation** transformation)
{
	_transformation = transformation;
	_origintransformation = irtkTransformation::New(transformation[0]);
}

inline void irtkCardiacSpatialCorrection::Debug(string message)
{
  if (_DebugFlag == true) cout << message << endl;
}

#endif
