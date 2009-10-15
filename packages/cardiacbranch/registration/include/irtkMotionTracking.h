/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMOTIONTRACKING_H

#define _IRTKMOTIONTRACKING_H

/**
 * Generic for motion tracking using image registration.
 *
 * This class implements a motion registration filter which takes an image
 * sequence and calculates the transformation which maps each time frame into
 * the coordinate system of the first time frame.
 *
 */

class irtkMotionTracking : public irtkRegistration
{

  /// Interface to input file stream
  friend istream& operator>> (istream&, irtkMotionTracking*);

  /// Interface to output file stream
  friend ostream& operator<< (ostream&, const irtkMotionTracking*);

protected:

  /// Image sequence
  irtkGreyImage *_image;

  /// Mask for valid image region over which to track
  irtkGreyImage *_mask;

  /// Transformation
  irtkTransformation *_transformation;

  /// Transformation in B-Spline format
  irtkBSplineFreeFormTransformation4D *_affd;

  /// Histogram
  irtkSimilarityMetric *_metric;

  /// Interpolator
  irtkInterpolateImageFunction *_interpolator;

  /// Optimizer
  irtkOptimizer *_optimizer;

  /// Blurring of image (in mm)
  double _Blurring[MAX_NO_RESOLUTIONS];

  /// Resolution of image (in mm)
  double _Resolution[MAX_NO_RESOLUTIONS][3];

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

  /// Optimization method for registration
  irtkOptimizationMethod _OptimizationMethod;

  /// Interpolation mode to use during resampling and registration
  irtkInterpolationMode _InterpolationMode;

  /// Control point spacing in the x-direction
  double _DX;

  /// Control point spacing in the y-direction
  double _DY;

  /// Control point spacing in the z-direction
  double _DZ;

  /// Convergence parameter for optimization based on change in similarity.
  double _Epsilon;

  /// Convergence parameter for optimization based on change in the transformation.
  double _Delta[MAX_NO_RESOLUTIONS];

  /// Speedup factor when calculating derivative
  double _SpeedupFactor;

  /// Debugging flag
  int    _DebugFlag;

  /// Source image domain which can be interpolated fast
  double _x1, _y1, _z1, _x2, _y2, _z2;

  /// Used as temporary memory for metric
  irtkSimilarityMetric *_tmpMetricA, *_tmpMetricB;

  /** Used as lookup table for transformed coordinates up to level n-1.  This
      lookup table needs to be calculated only once for each image resolution
      level. */
  float *_mffdLookupTable;

  /** Used as lookup table for transformed coordinates including level n. This
      lookup table needs to be calculated each time a control point has been
      modified. */
  float *_affdLookupTable;

  /** Used as lookup table for the contribution of each control point. This
      lookup table needs to be calculated only once. */
  float *_localLookupTable;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

  /// Initial set up for the registration at a multiresolution level
  virtual void Initialize(int);

  /// Final set up for the registration at a multiresolution level
  virtual void Finalize(int);

  /// Update lookup table
  virtual void UpdateLUT();

public:

  /// Constructor
  irtkMotionTracking();

  /// Destructor
  virtual ~irtkMotionTracking();

  /// Sets input for the registration filter
  virtual void SetInput (irtkGreyImage *, irtkGreyImage *);

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *) = 0;

  /// Runs the registration filter
  virtual void Run();

  /** Evaluates the similarity metric. This function evaluates the similarity
   *  metric of the registration by looping over the target image and
   *  interpolating the transformed source image while filling the joint
   *  histogram. This function returns the value of the similarity measure
   *  using Similarity().
   */
  virtual double Evaluate() = 0;

  /** Evaluates the registration. This function evaluates the registration by
   *  looping over the target image and interpolating the transformed source
   *  image while filling the joint histogram. This function returns the value
   *  of the similarity measure using Similarity(). This function uses the
   *  cached result of any previous call to Evaluate() and recalculates the
   *  similarity measure in the specified region of interest.
   */
  virtual double EvaluateDerivative(int, double);

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
  virtual SetMacro(OptimizationMethod, irtkOptimizationMethod);
  virtual GetMacro(OptimizationMethod, irtkOptimizationMethod);

};

inline void irtkMotionTracking::SetInput(irtkGreyImage *image, irtkGreyImage *mask)
{
  _image = image;
  _mask  = mask;
}

inline void irtkMotionTracking::Debug(string message)
{
  if (_DebugFlag == True) cout << message << endl;
}

#endif

