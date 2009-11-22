/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKOPTIMIZER_H

#define _IRTKOPTIMIZER_H

/// Forward declaration
class irtkRegistration;

/**
 * Generic class for optimization of voxel-based registration.
 *
 * This class implements an optimizer which maximises voxel similarity
 * measures as used in image registration. This is the abstract base class
 * which defines a common interface for arbitrary optimization filters. Each
 * derived class has to implement all abstract member functions.
 *
 */

class irtkOptimizer : public irtkObject
{

protected:

  /// Pointer to transformation
  irtkTransformation *_Transformation;

  /// Pointer to registration
  irtkRegistration *_Registration;

  /// Step size
  double _StepSize;

  /// Epsilon
  double _Epsilon;

  /// Storage for monitoring change in the transformation.
  double *transformationParams;

public:

  /// Constructor
  irtkOptimizer();

  /// Destructor
  virtual ~irtkOptimizer();

  /// Run the optimizer
  virtual double Run() = 0;

  /// Run the optimizer
  virtual void Run(double &, double &);

  /// Print name of the class
  virtual const char *NameOfClass() = 0;

  virtual void SetTransformation(irtkTransformation *);

  virtual GetMacro(Transformation, irtkTransformation *);

  virtual SetMacro(Registration, irtkRegistration *)

  virtual GetMacro(Registration, irtkRegistration *);

  virtual SetMacro(StepSize, double);

  virtual GetMacro(StepSize, double);

  virtual SetMacro(Epsilon, double);

  virtual GetMacro(Epsilon, double);

};

inline void irtkOptimizer::SetTransformation(irtkTransformation *transformation)
{
  int n;
  _Transformation = transformation;
  n = _Transformation->NumberOfDOFs();
  transformationParams = new double[n];
}

#include <irtkDownhillDescentOptimizer.h>
#include <irtkGradientDescentOptimizer.h>
#include <irtkSteepestGradientDescentOptimizer.h>
#include <irtkConjugateGradientDescentOptimizer.h>

#endif
