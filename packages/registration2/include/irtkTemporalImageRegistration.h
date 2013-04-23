/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKTEMPROALIMAGEREGISTRATION_H

#define _IRTKTEMPROALIMAGEREGISTRATION_H

/**
 * Generic for image registration based on voxel similarity measures.
 *
 * This class implements a registration filter which takes two input images
 * and calculates the transformation which maps the second image (denotes as
 * source image) into the coordinate system of the first image (denotes as
 * target image).  This is the abstract base class which defines a common
 * interface for arbitrary registration filters. Each derived class has to
 * implement all abstract member functions.
 *
 */

class irtkTemporalImageRegistration : public irtkRegistration2
{

  /// Interface to input file stream
  friend istream& operator>> (istream&, irtkTemporalImageRegistration*);

  /// Interface to output file stream
  friend ostream& operator<< (ostream&, const irtkTemporalImageRegistration*);

protected:

  /** Input image. This image is denoted as target image and its
   *  coordinate systems define the frames of reference for the registration.
   */
  irtkGenericImage<short> *_target;

  /** Input image sequence. This image is denoted as source images. The goal of
   *  the registration is to find the transformation which maps the source
   *  image sequence into the coordinate system of the target image.
   */
  irtkGenericImage<short> **_source;
  int _N_source;

  /** Current estimate of the source image transformed back into the target
   *  coordinate system. This is updated every time the Update function is
   *  called.
   */
  irtkGenericImage<double> *_transformedSource;

  /// Gradient of the original source
  irtkGenericImage<double> *_sourceGradient;

  /// Gradient of the transformed source
  irtkGenericImage<double> *_transformedSourceGradient;

  /// Gradient of the similarity metric
  irtkGenericImage<double> *_similarityGradient;

  /// Transformation
  irtkTransformation *_transformation;

  /// 2D histogram (this is not used for all similarity metrics)
  irtkHistogram_2D<double> *_histogram;

  /// Interpolator for source image
  irtkInterpolateImageFunction **_interpolator;

  /// Interpolator for source image gradient
  irtkInterpolateImageFunction **_interpolatorGradient;

  /// Blurring of target image (in mm)
  double _TargetBlurring[MAX_NO_RESOLUTIONS];

  /// Resolution of target image (in mm)
  double _TargetResolution[MAX_NO_RESOLUTIONS][3];

  /// Blurring of source image (in mm)
  double _SourceBlurring[MAX_NO_RESOLUTIONS];

  /// Resolution of source image (in mm)
  double _SourceResolution[MAX_NO_RESOLUTIONS][3];

  /// Minimum length of steps
  double _MinStep[MAX_NO_RESOLUTIONS];

  /// Maximum length of steps
  double _MaxStep[MAX_NO_RESOLUTIONS];

  /// Max. number of iterations per step size
  int    _NumberOfIterations[MAX_NO_RESOLUTIONS];

  /// Padding value of target image
  short  _TargetPadding;

  /// Padding value of source image
  short _SourcePadding;

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

  /// Debugging flag
  int    _DebugFlag;

  /// Current min and max voxel values
  int _target_min, _target_max;
  int _source_min, _source_max;
  int _maxDiff;

  /// Current level in the multi-resolution pyramid
  int _CurrentLevel;

  /// Current iteration during the optimization
  int _CurrentIteration;

  /// Source image domain which can be interpolated fast
  double _source_x1, _source_y1, _source_z1;
  double _source_x2, _source_y2, _source_z2;

  /// Dimension of subdivision (4 = 4D, 3 = 3D, 2 = 2D(x,y))
  int _SubdivisionDim[MAX_NO_RESOLUTIONS];

  /// Array containing the real (relative) time points of the targets
  double * _t_real;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

  /// Initial set up for the registration at a multiresolution level
  virtual void Initialize(int);

  /// Final set up for the registration at a multiresolution level
  virtual void Finalize(int);

  /// Update state of the registration based on current transformation estimate
  virtual void Update(bool);

  /// Update state of the registration based on current transformation estimate (source image)
  virtual void UpdateSource() = 0;

  /// Update state of the registration based on current transformation estimate (source image and source image gradient)
  virtual void UpdateSourceAndGradient() = 0;

  /// Evaluate similarity measure: SSD
  virtual double EvaluateSSD();

  /// Evaluate gradient of similarity measure: SSD
  virtual void EvaluateGradientSSD();

  /// Evaluate similarity measure: NMI
  virtual double EvaluateNMI();

  /// Evaluate gradient of similarity measure: NMI
  virtual void EvaluateGradientNMI();

public:

  /// Constructor
  irtkTemporalImageRegistration();

  /// Destructor
  virtual ~irtkTemporalImageRegistration();

  /// Sets input for the registration filter
  virtual void SetInput (irtkGreyImage *, irtkGreyImage **, int);

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *) = 0;

  /// Set array containing the real (relative) time points of the targets
  virtual void SetTime(double *);

  /// Runs the registration filter
  virtual void Run();

  /** Evaluates the similarity metric. This function evaluates the similarity
   *  metric of the registration by looping over the target image and
   *  interpolating the transformed source image. This function returns the
   *  value of the similarity measure.
   */
  virtual double Evaluate();

  /** Evaluates the gradient of the similarity metric. This function
   *  evaluates the gradient of the similarity metric of the registration
   *  by looping over the target image and interpolating the transformed
   *  source image while computing the gradient. The derivatives in the
   *  gradient are expressed with respect to voxel wise displacements, not
   *  with respect to the parameters of the transformation. Each subclass
   *  of this class must override this function in order to compute the
   *  gradient with respect to the transformation parameters.
   */
  virtual double EvaluateGradient(double *);

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
  virtual SetMacro(OptimizationMethod, irtkOptimizationMethod);
  virtual GetMacro(OptimizationMethod, irtkOptimizationMethod);

};

inline void irtkTemporalImageRegistration::SetInput(irtkGreyImage *target, irtkGreyImage **source, int N_source)
{
  _N_source = N_source;
  _source = new irtkGreyImage *[_N_source];
  for (int i=0; i<_N_source; i++)
    _source[i] = source[i];
  _target = target;
}

inline void irtkTemporalImageRegistration::SetTime(double * t)
{
  _t_real = t;
}

inline void irtkTemporalImageRegistration::Debug(string message)
{
  if (_DebugFlag == true) cout << message << endl;
}

#include <irtkImageTFFDRegistration.h>

#endif
