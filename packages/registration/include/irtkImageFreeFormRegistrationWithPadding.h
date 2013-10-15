/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id: irtkImageFreeFormRegistrationWithPadding.h 304 2011-04-10 21:46:30Z ws207 $
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date: 2011-04-10 22:46:30 +0100 (Sun, 10 Apr 2011) $
  Version   : $Revision: 304 $
  Changes   : $Author: ws207 $

=========================================================================*/

#ifndef _IRTKIMAGEFREEFORMREGISTRATIONWITHPADDING_H

#define _IRTKIMAGEFREEFORMREGISTRATIONWITHPADDING_H

// Definition of available registration modes
//typedef enum { RegisterX, RegisterY, RegisterXY, RegisterXYZ } irtkImageFreeFormRegistrationMode;

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

class irtkImageFreeFormRegistrationWithPadding : public irtkImageRegistrationWithPadding
{

#ifdef HAS_TBB

  friend class irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluate;
  friend class irtkMultiThreadedImageFreeFormRegistrationWithPaddingEvaluateGradient;

#endif

  /// Friend declaration of NR optimization routines
  friend float irtkFreeFormRegistration_Ptr2NRfunc (float *x);

  /// Friend declaration of NR optimization routines
  friend void  irtkFreeFormRegistration_Ptr2NRdfunc(float *x, float *dx);

protected:

  /// Pointer to the local transformation which is currently optimized
  irtkBSplineFreeFormTransformation *_affd;

  /// Pointer to the global transformation which is constant
  irtkMultiLevelFreeFormTransformation *_mffd;

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

  /// Smoothness parameter for non-rigid registration
  double _Lambda1;

  /// Volume preservation parameter for non-rigid registration
  double _Lambda2;

  /// Topology preservation parameter for non-rigid registration
  double _Lambda3;

  /// Control point spacing in the x-direction
  double _DX;

  /// Control point spacing in the y-direction
  double _DY;

  /// Control point spacing in the z-direction
  double _DZ;

  /// Subdivide FFD between resolution levels
  bool _Subdivision;

  /// Speedup factor when calculating derivative
  double _SpeedupFactor;

  /// Registration mode
  irtkImageFreeFormRegistrationMode _Mode;

  /// Multilevel mode
  bool _MFFDMode;

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

  /// Update lookup table
  virtual void UpdateLUT();

public:

  /// Constructor
  irtkImageFreeFormRegistrationWithPadding();

  /// Set output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Print some information
  virtual void Print();

  /// Guess parameters
  virtual void GuessParameter();

  /// Guess parameters for distortion correction
  virtual void GuessParameterDistortion(double resolution = 0, double CPS=0, double penalty = 0);

  /// Read single line of registration parameters
  virtual bool Read(char *, char *, int &);

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

  virtual SetMacro(MFFDMode, bool);
  virtual GetMacro(MFFDMode, bool);
};

inline void irtkImageFreeFormRegistrationWithPadding::SetOutput(irtkTransformation *transformation)
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistrationWithPadding::SetOutput");

  if (strcmp(transformation->NameOfClass(),
             "irtkMultiLevelFreeFormTransformation") != 0) {
    cerr << "irtkImageFreeFormRegistrationWithPadding::SetOutput: Transformation must be "
         << "irtkMultiLevelFreeFormTransformation" << endl;
    exit(0);
  }
  _transformation = transformation;
}

inline const char *irtkImageFreeFormRegistrationWithPadding::NameOfClass()
{
  return "irtkFreeFormRegistrationWithPadding";
}

inline void irtkImageFreeFormRegistrationWithPadding::Print()
{}

#endif
