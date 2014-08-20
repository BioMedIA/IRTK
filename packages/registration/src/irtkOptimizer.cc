/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#include <irtkRegistration.h>

#include <gsl/gsl_vector.h>

// Global variable used for NR optimization
irtkOptimizer *_optimizer;

irtkOptimizer::irtkOptimizer()
{
  _optimizer = this;

  // Initialuze variables
  _Epsilon  = 0.001;
  _StepSize = 0.1;
}

// Used for GSL optimization: Evaluates similarity measure
float irtkRegistrationEvaluate(float *x)
{
  int i;

  for (i = 0; i < _optimizer->GetTransformation()->NumberOfDOFs(); i++) {
    _optimizer->GetTransformation()->Put(i, x[i+1]);
  }
  return -_optimizer->GetRegistration()->Evaluate();
}

// Used for GSL optimization: Evaluates similarity measure
void irtkRegistrationEvaluateDerivative(float *x, float *dx)
{
  int i;

  for (i = 0; i < _optimizer->GetTransformation()->NumberOfDOFs(); i++) {
    _optimizer->GetTransformation()->Put(i, x[i+1]);
  }
  _optimizer->GetRegistration()->EvaluateGradient(_optimizer->GetStepSize(), dx+1);
}


/* The function f that we want to minimize */
double evaluate_f (const gsl_vector *v, void *params)
{
  size_t i;

  float *x = new float[v->size];
  for (i = 0; i < v->size; i++) {
     x[i] = gsl_vector_get(v, i);
  }

  double result = irtkRegistrationEvaluate(x-1);
  delete [] x;

  return result;
}

/* The gradient of f */
void evaluate_df (const gsl_vector *v, void *params, gsl_vector *df)
{
  size_t i;

  float *x = new float[v->size];
  float *dx = new float[v->size];

  for (i = 0; i < v->size; i++) {
     x[i] = gsl_vector_get(v, i);
  }

  irtkRegistrationEvaluateDerivative(x-1, dx-1);

  for (i = 0; i < v->size; i++) {
     gsl_vector_set(df, i, dx[i]);
  }

  delete [] x;
  delete [] dx;
}

/* Compute both f and df together */
void evaluate_fdf (const gsl_vector *x, void *params,
                   double *f, gsl_vector *df)
{
  *f = evaluate_f(x, params);
  evaluate_df(x, params, df);
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
