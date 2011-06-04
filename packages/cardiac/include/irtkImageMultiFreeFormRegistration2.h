/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGEMULTIFREEFORMREGISTRATION2_H

#define _IRTKIMAGEMULTIFREEFORMREGISTRATION2_H

/**
 * Filter for non-rigid registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the non-rigid registration
 * of two images given a segmentation.
 *
 */

class irtkImageMultiFreeFormRegistration2 : public irtkImageRegistration2
{

protected:
  /// Pointer to segmentation
  irtkGreyImage *_segmentation;

  /// Number of labels
  int _NumberOfLabels;

  /// multi transformations
  irtkTransformation **_mtransformation;
  /// Pointer to the local transformation which is currently optimized
  irtkBSplineFreeFormTransformation **_affd;

  /// Pointer to the global transformation which is constant
  irtkMultiLevelFreeFormTransformation **_mffd;

  /// Pointer to Adjugate jacobian matrix
  irtkMatrix **_adjugate;

  /// Pointer to jacobian determine
  double **_determine;

  /// Smoothness parameter for non-rigid registration
  double _Lambda1;

  /// Volume preservation parameter for non-rigid registration
  double _Lambda2;

  /// Topology preservation parameter for non-rigid registration
  double _Lambda3;

  /// Lambda coefficient for every segments
  double *_LambdaC;

  /// Control point spacing in the x-direction
  double _DX;

  /// Control point spacing in the y-direction
  double _DY;

  /// Control point spacing in the z-direction
  double _DZ;

  /// Subdivide FFD between resolution levels
  bool _Subdivision;

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

  /** Evaluates the volume preservation term. */
  virtual double VolumePreservationPenalty();

  /** Evaluates the volume preservation term. */
  virtual void VolumePreservationPenalty(int,int,double*);

  /// Update transformed source and source gradient
  virtual void Update();

  /** Evaluates the registration. This function evaluates the registration by
   *  looping over the target image and interpolating the transformed source
   *  image while filling the joint histogram. This function returns the value
   *  of the similarity measure using Similarity().
   */
  virtual double Evaluate();

  /// Evaluate the gradient of the similarity measure for the current transformation.
  virtual double EvaluateGradient(int, double *);

public:

  /// Constructor
  irtkImageMultiFreeFormRegistration2();

  /// Sets input for the registration filter
  virtual void SetInput (irtkGreyImage *, irtkGreyImage *, irtkGreyImage *);

  /// Set output for the registration filter
  virtual void SetOutput(irtkTransformation **);

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

  // Access parameters for control point space
  virtual SetMacro(DX, double);
  virtual GetMacro(DX, double);
  virtual SetMacro(DY, double);
  virtual GetMacro(DY, double);
  virtual SetMacro(DZ, double);
  virtual GetMacro(DZ, double);
  virtual SetMacro(LambdaC, double*);

  // Access parameters for registration mode
  virtual SetMacro(Mode, irtkImageFreeFormRegistrationMode);
  virtual GetMacro(Mode, irtkImageFreeFormRegistrationMode);

  virtual SetMacro(MFFDMode, bool);
  virtual GetMacro(MFFDMode, bool);
};

inline void irtkImageMultiFreeFormRegistration2::SetOutput(irtkTransformation *transformation)
{
    // Print debugging information
    this->Debug("irtkImageFreeFormRegistration::SetOutput");

        if (strcmp(transformation->NameOfClass(),
            "irtkMultiLevelFreeFormTransformation") != 0) {
                cerr << "irtkImageMultiFreeFormRegistration::SetOutput: Transformation must be "
                    << "irtkMultiLevelFreeFormTransformation" << endl;
                exit(0);
        }
        cerr << "irtkImageMultiFreeFormRegistration::SetOutput: Transformation must be "
            << "Multiple irtkMultiLevelFreeFormTransformation" << endl;
        _transformation = transformation;
}

inline const char *irtkImageMultiFreeFormRegistration2::NameOfClass()
{
  return "irtkImageMultiFreeFormRegistration";
}

inline void irtkImageMultiFreeFormRegistration2::Print()
{}

#endif
