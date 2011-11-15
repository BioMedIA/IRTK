/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKDISCONTINUOUSFREEFORMREGISTRATION_H

#define _IRTKDISCONTINUOUSFREEFORMREGISTRATION_H

/**
 * Filter for non-rigid registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the non-rigid registration
 * of two images. The basic algorithm is described in Rueckert et al., IEEE
 * Transactions on Medical Imaging, In press.
 *
 */

class irtkDiscontinuousFreeFormRegistration : public irtkImageRegistration2
{

protected:

  /// Pointer to the local transformation which is currently optimized
  irtkBSplineFreeFormTransformation *_affd;

  /// Pointer to the global transformation which is constant
  irtkMultiLevelFreeFormTransformation *_mffd;

  /// Pointer to the local transformation which is currently optimized
  irtkBSplineFreeFormTransformation *_bsfd;

  /// Pointer to the global transformation which is constant
  irtkMultiLevelFreeFormTransformation *_bssm;

  /// Gradient of the similarity metric
  irtkGenericImage<double> _gradientResidual;

  /// Pointer to adjugate Jacobian matrix
  irtkMatrix *_adjugate;

  /// Pointer to Jacobian determinant
  double *_determinant;

  /// current gradient
  double *_currentgradient;

  double _compressrate;

  double *_sparsityindex;

  int _currentffdlevel;

  /// Volume parameter for non-rigid registration
  double _Lambda1;

  /// Registration mode
  irtkImageFreeFormRegistrationMode _Mode;

  /// Number of sparse model levels
  int _NumberOfModels;

  /// Start levels
  int _LevelOfStart;

  /// Previous ffd number for conjugate
  int _previousFFD;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Initial set up for the registration
  virtual void Initialize(int);

  /// Final set up for the registration
  virtual void Finalize();

  /// Final set up for the registration
  virtual void Finalize(int);

  /// Update state of the registration based on current transformation estimate
  virtual void Update(bool);

  /// Update jocobian and determine matrixes
  virtual void UpdateVolume(bool);

  /** Evaluates the volume preservation term. */
  virtual double VolumePreservationPenalty();

  /** Evaluates the volume preservation term. */
  virtual void VolumePreservationPenaltyGradient();

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
  virtual double EvaluateGradient(double *);

  /// Project gradient to B Spline domian
  virtual void DomainDecomposition();

  /// Choose which level to use
  virtual void ChooseLevel();

public:

  /// Constructor
  irtkDiscontinuousFreeFormRegistration();

  ~irtkDiscontinuousFreeFormRegistration();

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

};

inline void irtkDiscontinuousFreeFormRegistration::SetOutput(irtkTransformation *transformation)
{
  // Print debugging information
  this->Debug("irtkImageFreeFormRegistration::SetOutput");

  if (strcmp(transformation->NameOfClass(),
             "irtkMultiLevelFreeFormTransformation") != 0) {
    cerr << "irtkImageFreeFormRegistration::SetOutput: Transformation must be "
         << "irtkMultiLevelFreeFormTransformation" << endl;
    exit(0);
  }
  _transformation = transformation;
}

inline const char *irtkDiscontinuousFreeFormRegistration::NameOfClass()
{
  return "irtkFreeFormRegistration";
}

inline void irtkDiscontinuousFreeFormRegistration::Print()
{}

#endif
