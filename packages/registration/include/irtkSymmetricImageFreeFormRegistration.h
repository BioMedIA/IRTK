/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSYMMETRICIMAGEFREEFORMREGISTRATION_H

#define _IRTKSYMMETRICIMAGEFREEFORMREGISTRATION_H

/**
 * Filter for non-rigid registration based on voxel similarity measures.
 *
 * This class implements a registration filter for the non-rigid registration
 * of two images. The basic algorithm is described in Rueckert et al., IEEE
 * Transactions on Medical Imaging, In press.
 *
 */

class irtkSymmetricImageFreeFormRegistration : public irtkSymmetricImageRegistration
{

  /// Friend declaration of NR optimization routines
  friend float irtkFreeFormRegistration_Ptr2NRfunc (float *x);

  /// Friend declaration of NR optimization routines
  friend void  irtkFreeFormRegistration_Ptr2NRdfunc(float *x, float *dx);

protected:

  /// Pointer to the local transformation which is currently optimized
  irtkBSplineFreeFormTransformation *_affd1;
  irtkBSplineFreeFormTransformation *_affd2;

  /// Pointer to the global transformation which is constant
  irtkMultiLevelFreeFormTransformation *_mffd1;
  irtkMultiLevelFreeFormTransformation *_mffd2;

  /// Used as temporary memory for metric
  irtkSimilarityMetric *_tmp1MetricA, *_tmp1MetricB;
  irtkSimilarityMetric *_tmp2MetricA, *_tmp2MetricB;

  /** Used as lookup table for transformed coordinates up to level n-1.  This
      lookup table needs to be calculated only once for each image resolution
      level. */
  float *_mffdLookupTable1;
  float *_mffdLookupTable2;

  /** Used as lookup table for transformed coordinates including level n. This
      lookup table needs to be calculated each time a control point has been
      modified. */
  float *_affdLookupTable1;
  float *_affdLookupTable2;

  /** Used as lookup table for the contribution of each control point. This
      lookup table needs to be calculated only once. */
  float *_localLookupTable;

  /// Smoothness parameter for non-rigid registration
  double _Lambda1;

  /// Volume preservation parameter for non-rigid registration
  double _Lambda2;

  /// Topology preservation parameter for non-rigid registration
  double _Lambda3;

  /// Inverse consistency parameter for non-rigid registration
  double _Lambda4;

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

  /** Evaluates the inverse consistency term. */
  virtual double InverseConsistencyPenalty();

  /** Evaluates the topology preservation term. */
  virtual double InverseConsistencyPenalty(int);

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
  virtual double EvaluateDerivative1(int, double);
  virtual double EvaluateDerivative2(int, double);

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
  irtkSymmetricImageFreeFormRegistration();

  /// Set output for the registration filter
  virtual void SetOutput(irtkTransformation *, irtkTransformation *);

  /// Returns the name of the class
  virtual const char *NameOfClass();

  /// Print some information
  virtual void Print();

  /// Guess parameters
  virtual void GuessParameter();

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
};

inline void irtkSymmetricImageFreeFormRegistration::SetOutput(irtkTransformation *transformation1,
    irtkTransformation *transformation2)
{
  // Print debugging information
  this->Debug("irtkSymmetricImageFreeFormRegistration::SetOutput");

  if (strcmp(transformation1->NameOfClass(),
             "irtkMultiLevelFreeFormTransformation") != 0) {
    cerr << "irtkSymmetricImageFreeFormRegistration::SetOutput: Transformation must be "
         << "irtkMultiLevelFreeFormTransformation" << endl;
    exit(0);
  }
  if (strcmp(transformation2->NameOfClass(),
             "irtkMultiLevelFreeFormTransformation") != 0) {
    cerr << "irtkSymmetricImageFreeFormRegistration::SetOutput: Transformation must be "
         << "irtkMultiLevelFreeFormTransformation" << endl;
    exit(0);
  }
  _transformation1 = transformation1;
  _transformation2 = transformation2;
}

inline const char *irtkSymmetricImageFreeFormRegistration::NameOfClass()
{
  return "irtkFreeFormRegistration";
}

inline void irtkSymmetricImageFreeFormRegistration::Print()
{}

#endif
