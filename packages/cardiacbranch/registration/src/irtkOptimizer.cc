/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

// Global variable used for NR optimization
irtkOptimizer *_optimizer;

irtkOptimizer::irtkOptimizer()
{
  _optimizer = this;

  // Initialuze variables
  _Epsilon  = 0.001;
  _StepSize = 0.1;
}

// Used for NR optimization: Evaluates similarity measure
float irtkRegistrationEvaluate(float *x)
{
  int i;

  for (i = 0; i < _optimizer->GetTransformation()->NumberOfDOFs(); i++) {
    _optimizer->GetTransformation()->Put(i, x[i+1]);
  }
  return -_optimizer->GetRegistration()->Evaluate();
}

// Used for NR optimization: Evaluates similarity measure
void irtkRegistrationEvaluateDerivative(float *x, float *dx)
{
  int i;

  for (i = 0; i < _optimizer->GetTransformation()->NumberOfDOFs(); i++) {
    _optimizer->GetTransformation()->Put(i, x[i+1]);
  }
  _optimizer->GetRegistration()->EvaluateGradient(_optimizer->GetStepSize(), dx+1);
}

irtkOptimizer::~irtkOptimizer()
{
  delete [] transformationParams;
}

void irtkOptimizer::Run(double &epsilon, double &maxChange)
{
  int i, n;
  double diff;
  n = _Transformation->NumberOfDOFs();

  // Store the current state of the transformation.
  for (i = 0; i < n; ++i) {
    transformationParams[i] = _Transformation->Get(i);
  }

  // Epsilon is the change in similarity over an iteration of the
  // optimiser.
  epsilon = this->Run();

  // MaxChange is the maximum change over the transformation parameters.
  maxChange = 0;
  for (i = 0; i < n; ++i) {
    diff = fabs(_Transformation->Get(i) - transformationParams[i]);

    if (maxChange < diff) {
      maxChange = diff;
    }
  }
}
