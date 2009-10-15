#ifndef _IRTKCARDIACIMAGEFREEFORMREGISTRATION_H

#define _IRTKCARDIACIMAGEFREEFORMREGISTRATION_H

/// MaxWeight for the Untagged and Tagged Images
#define MaxWeight 3

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

class irtkCardiacImageFreeFormRegistration : public irtkImageRegistration
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
  /** First input taggedimage. This image is denoted as taggedtarget image and its
   *  coordinate system defines the frame of reference for the registration.
   */
  irtkGreyImage *_ttarget;

  /** Second input taggedimage. This image is denoted as taggedsource image. The goal of
   *  the registration is to find the transformation which maps the source
   *  image into the coordinate system of the target image.
   */
  irtkGreyImage *_tsource;

   /// tagged source's Interpolator
  irtkInterpolateImageFunction *_tinterpolator;

  /// Pointer to the local transformation which is currently optimized
  irtkBSplineFreeFormTransformation *_affd;

  /// Pointer to the global transformation which is constant
  irtkMultiLevelFreeFormTransformation *_mffd;

  /// Used as temporary memory for metric
  irtkSimilarityMetric *_tmpMetricA, *_tmpMetricB;

  /// Used as temporary memory for metric
  irtkSimilarityMetric *_ttmpMetricA, *_ttmpMetricB;

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

  /// Smoothness parameter for non-rigid registration
  double _Lambda1;

  /// Volume preservation parameter for non-rigid registration
  double _Lambda2;

  /// Topology preservation parameter for non-rigid registration
  double _Lambda3;

  /// Weight parameter for tagged images
  double _Weight;

  /// Weight parameter for tagged images
  double _Tweight;

  /// Threshold to identify the heart wall
  double _Thresholdmax;

  /// Threshold to identify the heart wall
  double _Thresholdmin;

  /// Control point spacing in the x-direction
  double _DX;

  /// Control point spacing in the y-direction
  double _DY;

  /// Control point spacing in the z-direction
  double _DZ;

  /// Subdivide FFD between resolution levels
  Bool _Subdivision;

  /// Speedup factor when calculating derivative
  double _SpeedupFactor;
  
  /// Registration mode
  irtkImageFreeFormRegistrationMode _Mode;

  double _tsource_x1, _tsource_y1, _tsource_z1;
  double _tsource_x2, _tsource_y2, _tsource_z2; 
  /// Initial set up for the registration
  virtual void Initialize();

  /// Initial set up for the registration
  virtual void Initialize(int);

  /// Final set up for the registration
  virtual void Finalize();

  /// Final set up for the registration
  virtual void Finalize(int);

  /** Evaluates the smoothness preservation term. */
  virtual double SmoothnessPenalty();

  /** Evaluates the smoothness term. */
  virtual double SmoothnessPenalty(int);

  /** Evaluates the volume preservation term. */
  virtual double VolumePreservationPenalty();

  /** Evaluates the volume preservation term. */
  virtual double VolumePreservationPenalty(int);

  /** Evaluates the topology preservation term. */
  virtual double TopologyPreservationPenalty();

  /** Evaluates the topology preservation term. */
  virtual double TopologyPreservationPenalty(int);

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

  virtual void CorrespondTaggedImages(irtkGreyImage *, irtkGreyImage *);

  virtual void EvaluateWeight(irtkGreyImage *, irtkGreyImage *);

  virtual int WeightFunction (double value);

  /// Update lookup table
  virtual void UpdateLUT();

public:

  /// Constructor
  irtkCardiacImageFreeFormRegistration();

  /// Set output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Print some information
  virtual void Print();

  /// Guess parameters
  virtual void GuessParameter();

  /// Sets input for the registration filter
  virtual void SetInput (irtkGreyImage *, irtkGreyImage *, irtkGreyImage *, irtkGreyImage *);

  /// Read single line of registration parameters
  virtual Bool Read(char *, char *, int &);

  /// Write registration parameters to file
  virtual void Write(ostream &);

  // Access parameters for control point space
  virtual SetMacro(DX, double);
  virtual GetMacro(DX, double);
  virtual SetMacro(DY, double);
  virtual GetMacro(DY, double);
  virtual SetMacro(DZ, double);
  virtual GetMacro(DZ, double);
  virtual SetMacro(SpeedupFactor, double);
  virtual GetMacro(SpeedupFactor, double);

  // Access parameters for registration mode
  virtual SetMacro(Mode, irtkImageFreeFormRegistrationMode);
  virtual GetMacro(Mode, irtkImageFreeFormRegistrationMode);
};

inline void irtkCardiacImageFreeFormRegistration::SetOutput(irtkTransformation *transformation)
{
  // Print debugging information
  this->Debug("irtkCardiacImageFreeFormRegistration::SetOutput");

  if (strcmp(transformation->NameOfClass(),
             "irtkMultiLevelFreeFormTransformation") != 0) {
    cerr << "irtkCardiacImageFreeFormRegistration::SetOutput: Transformation must be "
         << "irtkMultiLevelFreeFormTransformation" << endl;
    exit(0);
  }
  _transformation = transformation;
}

inline const char *irtkCardiacImageFreeFormRegistration::NameOfClass()
{
  return "irtkCardiacImageFreeFormRegistration";
}

inline void irtkCardiacImageFreeFormRegistration::Print()
{}

#endif
