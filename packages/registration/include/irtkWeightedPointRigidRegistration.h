
#ifndef _IRTKWEIGHTEDPOINTRIGIDREGISTRATION_H

#define _IRTKWEIGHTEDPOINTRIGIDREGISTRATION_H

#include <irtkImage.h>

#include <irtkTransformation.h>

/**
 * Filter for point-based registration.
 *
 * This class implements a registration filter for point-based registration
 * of two sets of points. As an output, it returns the rigid transformation
 * coefficients in a irtkRigidTransformation class.
 *
 * every point pair has a weight associated with it
 *
*/

class irtkWeightedPointRigidRegistration : public irtkPointRigidRegistration
{

protected:

  /// point set weights
  double * _weights;

  /// Initial set up for the registration
  virtual void Initialize();

  // Final clean up
  virtual void Finalize();

  /// Optimize registration using closed form solution
  virtual void ClosedFormOptimizer();

public:

  /// Constructor
  irtkWeightedPointRigidRegistration();

  /// Destructor
  virtual ~irtkWeightedPointRigidRegistration();

  /// Sets output for the registration filter
  virtual void SetPointWeights(double *);

  /// Returns the name of the class
  virtual const char *NameOfClass();

};

#endif
