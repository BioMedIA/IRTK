/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCARDIAC3DIMAGEFREEFORMREGISTRATION_H

#define _IRTKCARDIAC3DIMAGEFREEFORMREGISTRATION_H

#ifdef HAS_TBB

class irtkMultiThreadedImageFreeFormRegistrationEvaluate;
class irtkMultiThreadedImageFreeFormRegistrationEvaluateGradient;

#endif

/**
 * Filter for non-rigid registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the non-rigid registration
 * of two images. The basic algorithm is described in Rueckert et al., IEEE
 * Transactions on Medical Imaging, In press.
 *
 */

class irtkCardiac3DImageFreeFormRegistration : public irtkMultipleImageFreeFormRegistration
{

#ifdef HAS_TBB

  friend class irtkMultiThreadedImageFreeFormRegistrationEvaluate;
  friend class irtkMultiThreadedImageFreeFormRegistrationEvaluateGradient;

#endif

  /// Friend declaration of NR optimization routines
  friend float irtkFreeFormRegistration_Ptr2NRfunc (float *x);

  /// Friend declaration of NR optimization routines
  friend void  irtkFreeFormRegistration_Ptr2NRdfunc(float *x, float *dx);

protected:

  /** First set of input image. This image is denoted as target image and its
   *  coordinate system defines the frame of reference for the registration.
   */
  irtkGreyImage **_utarget;

  /** Second input image. This image is denoted as source image. The goal of
   *  the registration is to find the transformation which maps the source
   *  image into the coordinate system of the target image.
   */
  irtkGreyImage **_usource;

  /// segmentation image
  irtkGreyImage *_threshold;

  /// adaptive incompressible constraint
  irtkRealImage *_omega;

  /// omega for control points
  double *_comega;

  /// myocardium probability
  irtkRealImage *_myoprob;

  /// Number of images (must be equal for source and target images)
  int _numberOfuImages;

  /// Histogram
  irtkSimilarityMetric **_umetric;

  /// Interpolator
  irtkInterpolateImageFunction **_uinterpolator;

   /// Blurring of target image (in mm)
  double _uTargetBlurring[MAX_NO_RESOLUTIONS];

  /// Resolution of target image (in mm)
  double _uTargetResolution[MAX_NO_RESOLUTIONS][3];

  /// Blurring of source image (in mm)
  double _uSourceBlurring[MAX_NO_RESOLUTIONS];

  /// Resolution of source image (in mm)
  double _uSourceResolution[MAX_NO_RESOLUTIONS][3];

  /// Source image domain which can be interpolated fast
  double *_usource_x1, *_usource_y1, *_usource_z1;
  double *_usource_x2, *_usource_y2, *_usource_z2;

  /// Used as temporary memory for metric
  irtkSimilarityMetric **_utmpMetricA, **_utmpMetricB;

  /** Used as lookup table for transformed coordinates up to level n-1.  This
      lookup table needs to be calculated only once for each image resolution
      level. */
  float **_umffdLookupTable;

  /** Used as lookup table for transformed coordinates including level n. This
      lookup table needs to be calculated each time a control point has been
      modified. */
  float **_uaffdLookupTable;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Initial set up for the registration
  virtual void Initialize(int);

  /// Final set up for the registration
  virtual void Finalize();

  /// Final set up for the registration
  virtual void Finalize(int);

  /** Evaluates the volume preservation term. */
  virtual double VolumePreservationPenalty();

  /** Evaluates the volume preservation term. */
  virtual double VolumePreservationPenalty(int);

  /** Evaluates the smoothness preservation term. */
  virtual double SmoothnessPenalty();

  /** Evaluates the smoothness term. */
  virtual double SmoothnessPenalty(int);

  /** Evaluates the omega for control point. */
  virtual void EvaluateOmega();

  /** Evaluates the omega for control point. */
  virtual void EvaluateOmega(int);

  /** Evaluates the registration. This function evaluates the registration by
   *  looping over the target image and interpolating the transformed source
   *  image while filling the joint histogram. This function returns the value
   *  of the similarity measure using Similarity().
   */
  virtual double Evaluate();

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

  /// Update lookup table
  virtual void UpdateLUT();
  /// Update adaptive jacobian constraint
  virtual void UpdateOmega();
  /// Evaluate weight functions
  virtual void EvaluateWeight(irtkRealImage *, irtkGreyImage *, irtkGreyImage *);
  /// Evaluate Myocardium Probability part one
  virtual void EvaluateMyoProb1(irtkRealImage *weight, irtkGreyImage *threshold, irtkGaussian &gaussian, double &denom);
  /// Evaluate Myocardium Probability part two
  virtual void EvaluateMyoProb2(irtkRealImage *weight, irtkGreyImage *threshold, irtkGaussian &gaussian, double &denom);
  /// Dialate the threshold
  virtual void ExtendThreshold(irtkGreyImage *, int);
  /// Evaluate the weight on a given voxel voxel
  virtual double WeightFunction (double edge, double edgemax, double threshold);

public:

  /// Constructor
  irtkCardiac3DImageFreeFormRegistration();

  /// Destructor
  virtual ~irtkCardiac3DImageFreeFormRegistration();

  /// Sets input for the registration filter
  virtual void SetInput (irtkGreyImage **, irtkGreyImage **, int, irtkGreyImage **, irtkGreyImage **, int);

  /// Sets input for the registration filter
  virtual void SetThreshold (irtkGreyImage *);

  /// Guess parameters
  virtual void GuessParameter();

    /// Runs the registration filter
  virtual void Run();

  /// Parse parameter line
  virtual bool Read(char *, char *, int &);

  /// Write parameters to stream
  virtual void Write(ostream &);

  /// Returns the name of the class
  virtual const char *NameOfClass();


};

inline const char *irtkCardiac3DImageFreeFormRegistration::NameOfClass()
{
  return "irtkCardiac3DImageFreeFormRegistration";
}

inline void irtkCardiac3DImageFreeFormRegistration::SetInput(irtkGreyImage **target, irtkGreyImage **source, int n1, 
															 irtkGreyImage **utarget, irtkGreyImage **usource, int n2)
{
	this->irtkMultipleImageFreeFormRegistration::SetInput(target,source,n1);	
	
	int i;
	_numberOfuImages = n2;
	_utarget = new irtkGreyImage*[_numberOfuImages];
	_usource = new irtkGreyImage*[_numberOfuImages];
	for (i = 0; i < _numberOfuImages; i++) {
		_utarget[i] = utarget[i];
		_usource[i] = usource[i];
	}
}

inline void irtkCardiac3DImageFreeFormRegistration::SetThreshold(irtkGreyImage *threshold)
{
	_threshold = threshold;
}

#endif

