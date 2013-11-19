/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGETSFFDREGISTRATION_H

#define _IRTKIMAGETSFFDREGISTRATION_H

/**
 * Filter for non-rigid registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the non-rigid registration
 * of two images. The basic algorithm is described in Rueckert et al., IEEE
 * Transactions on Medical Imaging, In press.
 *
 */

class irtkImageTSFFDRegistration : public irtkTemporalImageRegistration
{

protected:

  /// Pointer to the local transformation which is currently optimized
  irtkBSplineFreeFormTransformationPeriodic *_affd;

  /// Pointer to the global transformation which is constant
  irtkMultiLevelFreeFormTransformation *_mffd;

  /// variable for periodic transformation
  int _periodic;

  /// Pointer to lattice coordinates for every voxel
  double **_latticeCoordLUT;

  /// current gradient
  double *_currentgradient;

  double *_g;
  double *_h;
  double *_tmp;
  bool *_mask;

  int _NumberOfDofs;

  /// Number of sparse model levels
  int _NumberOfModels;

  /// Largest ffd grid spacing
  double _LargestSpacing;

  /// Finest ffd grid spacing
  double _FinestSpacing;

  /// Finest time grid spacing * interval
  double _FinestTimeSpacing;

  /// Marxium similarity
  double _MaxSimilarity;

  /// Current similarity
  double _CurrentSimilarity;

  /// Smoothness parameter for non-rigid registration
  double _Lambda1;

  /// temp for sparsity constrain
  double _Lambda2;

  /// Sparsity constraint lambda
  double _Lambda3;

  /// Sparsity constraint lambda
  bool _Lambda3off;

  /// Registration mode
  irtkImageFreeFormRegistrationMode _Mode;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Initial set up for the registration
  virtual void Initialize(int);

  /// Initialize transformation for registration level
  virtual void InitializeTransformation();

  /// Initialize transformation for registration level
  virtual void InitializeTransformation(int);

  /// Final set up for the registration
  virtual void Finalize();

  /// Final set up for the registration
  virtual void Finalize(int);

  /// Update state of the registration based on current transformation estimate
  virtual void Update(bool);

  /// Update state of the registration based on current transformation estimate (source image)
  virtual void UpdateSource();

  /// Update state of the registration based on current transformation estimate (source image and source image gradient)
  virtual void UpdateSourceAndGradient();

    /** Evaluates the smoothness preservation term. */
  virtual double SmoothnessPenalty();

  /** Evaluates the gradient of the smoothness term. */
  virtual void SmoothnessPenaltyGradient();

  /** Evaluates the sparse preservation term. */
  virtual double SparsePenalty();

  /** Evaluates the gradient of the sparse term. */
  virtual void SparsePenaltyGradient();

  /** Initialize the coordinate local cache */
  virtual void InitializeCoordLut();

  /** Initialize sparsity parameter with zero convergence */
  virtual void InitilizeSparsityParmeter();

  /** Evaluates the registration. This function evaluates the registration by
   *  looping over the target image and interpolating the transformed source
   *  image while filling the joint histogram. This function returns the value
   *  of the similarity measure using Similarity().
   */
  virtual double Evaluate();

  /// Evaluate the gradient of the similarity measure for the current transformation for 2D images
  virtual void EvaluateGradient2D();

  /// Evaluate the gradient of the similarity measure for the current transformation for 3D images
  virtual void EvaluateGradient3D();

  /// Evaluate the gradient of the similarity measure for the current transformation.
  virtual double EvaluateGradient();

  /// Normalization of similarity gradient
  virtual void NormalizeGradient();

public:

  /// Constructor
  irtkImageTSFFDRegistration();

  /// Set output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Runs the registration filter
  virtual void Run();

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Print some information
  virtual void Print();

  /// Guess parameters
  virtual void GuessParameter();

  /// Read single line of registration parameters
  virtual bool Read(char *, char *, int &);

  /// Write registration parameters to file
  virtual void Write(ostream &);

  // Access parameters for registration mode
  virtual SetMacro(Mode, irtkImageFreeFormRegistrationMode);
  virtual GetMacro(Mode, irtkImageFreeFormRegistrationMode);

  virtual SetMacro(Lambda1, double);
  virtual GetMacro(Lambda1, double);
  virtual SetMacro(Lambda3, double);
  virtual GetMacro(Lambda3, double);
};

inline void irtkImageTSFFDRegistration::SetOutput(irtkTransformation *transformation)
{
  // Print debugging information
  this->Debug("irtkImageTSFFDRegistration::SetOutput");

  if (strcmp(transformation->NameOfClass(),
             "irtkMultiLevelFreeFormTransformation") != 0) {
    cerr << "irtkImageTSFFDRegistration::SetOutput: Transformation must be "
         << "irtkMultiLevelFreeFormTransformation" << endl;
    exit(0);
  }
  _transformation = transformation;
}

inline const char *irtkImageTSFFDRegistration::NameOfClass()
{
  return "irtkImageTSFFDRegistration";
}

inline void irtkImageTSFFDRegistration::Print()
{}

#endif
