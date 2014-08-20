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

#include <gsl/gsl_multimin.h>

#define MAXITS 200

/* The function f that we want to minimize */
extern double evaluate_f (const gsl_vector *v, void *params);

/* The gradient of f */
extern void evaluate_df (const gsl_vector *v, void *params, gsl_vector *df);

/* Compute both f and df together */
extern void evaluate_fdf (const gsl_vector *x, void *params,
                   double *f, gsl_vector *df);

double irtkSteepestGradientDescentOptimizer::Run()
{
  int i, n, status;
  size_t iter = 0;
  float new_similarity, old_similarity;

  // Number of variables we have to optimize
  n = _Transformation->NumberOfDOFs();

  // Current similarity
  old_similarity = _Registration->Evaluate();
  cout << old_similarity << endl;

  // Convert some stuff to GSL
  gsl_vector *x = gsl_vector_alloc(n);
  for (i = 0; i < n; i++) {
    gsl_vector_set(x, i, _Transformation->Get(i));
  }

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_multimin_function_fdf func;
  func.n = n;
  func.f = &evaluate_f;
  func.df = &evaluate_df;
  func.fdf = &evaluate_fdf;

  T = gsl_multimin_fdfminimizer_steepest_descent;
  s = gsl_multimin_fdfminimizer_alloc(T, n);

  gsl_multimin_fdfminimizer_set(s, &func, x, _StepSize, _Epsilon);

  // Call GSL steepest descent routine
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status)
      break;

    status = gsl_multimin_test_gradient (s->gradient, _Epsilon);
  } while (status == GSL_CONTINUE && iter < MAXITS);

  new_similarity = s->f;

  // Convert some stuff back from GSL
  for (i = 0; i < n; i++) {
    _Transformation->Put(i, gsl_vector_get(s->x, i));
  }
  cout << -new_similarity << " " << endl;

  // Delete GSL memory
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  // Return
  if (-new_similarity > old_similarity) {
    return -new_similarity - old_similarity;
  } else {
    return 0;
  }
}
