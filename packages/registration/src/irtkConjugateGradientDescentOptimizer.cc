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

#include <gsl/gsl_multimin.h>

#define MAXITS 200

// Used for NR optimization: Evaluates similarity measure
extern float irtkRegistrationEvaluate(float *x);

// Used for NR optimization: Evaluates similarity measure
extern void  irtkRegistrationEvaluateDerivative(float *x, float *dx);

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

double irtkConjugateGradientDescentOptimizer::Run()
{
  int i, n, status;
  size_t iter = 0;
  float new_similarity, old_similarity;

  // Number of variables we have to optimize
  n = _Transformation->NumberOfDOFs();

  // Current similarity
  old_similarity = _Registration->Evaluate();

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

  T = gsl_multimin_fdfminimizer_conjugate_pr;
  s = gsl_multimin_fdfminimizer_alloc(T, n);

  gsl_multimin_fdfminimizer_set(s, &func, x, _StepSize, _Epsilon);

  // Call GSL Polak-Ribiere routine for conjugate gradient
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);

    if (status)
      break;

    status = gsl_multimin_test_gradient (s->gradient, _Epsilon);
  } while (status == GSL_CONTINUE && iter < MAXITS); 
  
  new_similarity = s->f;

  // Convert some stuff back from NR
  for (i = 0; i < n; i++) {
    _Transformation->Put(i, gsl_vector_get(s->x, i));
  }

  // Delete GSL memory
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  // Return
  if (-new_similarity > old_similarity) {
    cout << new_similarity << endl;
    return -new_similarity - old_similarity;
  } else {
    return 0;
  }
}
