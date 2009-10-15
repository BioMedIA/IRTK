/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKPOINTREGISTRATION_H

#define _IRTKPOINTREGISTRATION_H

#include <irtkRegistration.h>

#ifdef HAS_VTK
#include <vtkPolyData.h>
#endif

/**
 * Filter for point-based registration.
 *
 * This class implements a registration filter for point-based registration
 * of two sets of points.
 *
*/

class irtkPointRegistration : public irtkRegistration
{

protected:

  /// Input
  irtkPointSet *_target;

  /// Input
  irtkPointSet *_source;

  /// Output
  irtkTransformation *_transformation;

  /// Optimizer
  irtkOptimizer *_optimizer;

  /// Optimization method for registration
  irtkOptimizationMethod _OptimizationMethod;

  /// Initial set up for the registration
  virtual void Initialize();

  /// Final set up for the registration
  virtual void Finalize();

public:

  /// Constructor
  irtkPointRegistration();

  /// Destructor
  virtual ~irtkPointRegistration();

  /// Sets input for the registration filter
  virtual void SetInput (irtkPointSet *, irtkPointSet *);

  /// Sets output for the registration filter
  virtual void SetOutput(irtkTransformation *);

  /// Evaluates the similarity metric
  virtual double Evaluate();

  /// Evaluates the gradient of the similarity metric
  virtual double EvaluateGradient(float, float *);

  /// Run the filter
  virtual void Run() = 0;

  /// Returns the name of the class
  virtual const char *NameOfClass();

  // Access optimizer parameters
  virtual SetMacro(OptimizationMethod, irtkOptimizationMethod);
  virtual GetMacro(OptimizationMethod, irtkOptimizationMethod);
};

#include <irtkPointRigidRegistration.h>
#include <irtkPointAffineRegistration.h>
#include <irtkPointFreeFormRegistration.h>

#endif
