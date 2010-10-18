/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKIMAGEFLUIDREGISTRATION_H

#define _IRTKIMAGEFLUIDREGISTRATION_H

/**
 * Filter for non-rigid registration based on voxel similarity measures.
 *
 */


class irtkImageFluidRegistration : public irtkImageRegistration
{

  /// Friend declaration of NR optimization routines
  friend float irtkFluidRegistration_Ptr2NRfunc (float *x);

  /// Friend declaration of NR optimization routines
  friend void  irtkFluidRegistration_Ptr2NRdfunc(float *x, float *dx);

protected:

  /// Pointer to the local transformation which is currently optimized
  irtkBSplineFreeFormTransformation *_affd;

  /// Pointer to the global transformation which is constant
  irtkFluidFreeFormTransformation *_mffd;

  /// Control point spacing in the x-direction
  double _DX;

  /// Control point spacing in the y-direction
  double _DY;

  /// Control point spacing in the z-direction
  double _DZ;

  /// Number of time steps for fluid registration
  int _NumberOfTimeSteps;

  /// Speedup factor when calculating derivative
  double _SpeedupFactor;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Initial set up for the registration
  virtual void Initialize(int);

  /// Final set up for the registration
  virtual void Finalize();

  /// Final set up for the registration
  virtual void Finalize(int);

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

  virtual void UpdateLUT();

public:

  /// Constructor
  irtkImageFluidRegistration();

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

};

inline void irtkImageFluidRegistration::SetOutput(irtkTransformation *transformation)
{
  // Print debugging information
  this->Debug("irtkImageFluidRegistration::SetOutput");

  if (strcmp(transformation->NameOfClass(),
             "irtkFluidFreeFormTransformation") != 0) {
    cerr << "irtkImageFluidRegistration::SetOutput: Transformation must be "
         << "irtkFluidFreeFormTransformation" << endl;
    exit(0);
  }
  _transformation = transformation;
}

inline const char *irtkImageFluidRegistration::NameOfClass()
{
  return "irtkImageFluidRegistration";
}

inline void irtkImageFluidRegistration::Print()
{
}

#endif
