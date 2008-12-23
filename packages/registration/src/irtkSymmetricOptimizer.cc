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
irtkSymmetricOptimizer *_symmetricOptimizer;

irtkSymmetricOptimizer::irtkSymmetricOptimizer()
{
  _symmetricOptimizer = this;

  // Initialuze variables
  _Epsilon  = 0.001;
  _StepSize = 0.1;
}

// Used for NR optimization: Evaluates similarity measure
float irtkSymmetricRegistrationEvaluate(float *x)
{
  int i;

  for (i = 0; i < _symmetricOptimizer->GetTransformation1()->NumberOfDOFs(); i++) {
    _symmetricOptimizer->GetTransformation1()->Put(i, x[i+1]);
  }
  for (i = 0; i < _symmetricOptimizer->GetTransformation2()->NumberOfDOFs(); i++) {
    _symmetricOptimizer->GetTransformation2()->Put(i, x[i+_symmetricOptimizer->GetTransformation1()->NumberOfDOFs()+1]);
  }
  return -_symmetricOptimizer->GetRegistration()->Evaluate();
}

// Used for NR optimization: Evaluates similarity measure
void irtkSymmetricRegistrationEvaluateDerivative(float *x, float *dx)
{
  int i;

  for (i = 0; i < _symmetricOptimizer->GetTransformation1()->NumberOfDOFs(); i++) {
    _symmetricOptimizer->GetTransformation1()->Put(i, x[i+1]);
  }
  for (i = 0; i < _symmetricOptimizer->GetTransformation2()->NumberOfDOFs(); i++) {
    _symmetricOptimizer->GetTransformation2()->Put(i, x[i+_symmetricOptimizer->GetTransformation1()->NumberOfDOFs()+1]);
  }
  _symmetricOptimizer->GetRegistration()->EvaluateGradient(_symmetricOptimizer->GetStepSize(), dx+1);
}

irtkSymmetricOptimizer::~irtkSymmetricOptimizer()
{
  delete [] transformationParams;
}

void irtkSymmetricOptimizer::Run(double &epsilon, double &delta)
{
  int i;
  double diff;

  // Store the current state of the transformation.
  for (i = 0; i < _Transformation1->NumberOfDOFs(); ++i) {
    transformationParams[i] = _Transformation1->Get(i);
  }
  for (i = 0; i < _Transformation2->NumberOfDOFs(); ++i) {
    transformationParams[i+_Transformation1->NumberOfDOFs()] = _Transformation2->Get(i);
  }

  // Epsilon is the change in similarity over an iteration of the
  // optimiser.
  epsilon = this->Run();

  // Delta is the maximum change over the transformation parameters.
  delta = 0;
  for (i = 0; i < _Transformation1->NumberOfDOFs(); ++i) {
    diff = fabs(_Transformation1->Get(i) - transformationParams[i]);

    if (delta < diff) {
      delta = diff;
    }
  }
  for (i = 0; i < _Transformation2->NumberOfDOFs(); ++i) {
    diff = fabs(_Transformation2->Get(i) - transformationParams[i+_Transformation1->NumberOfDOFs()]);

    if (delta < diff) {
      delta = diff;
    }
  }
}
