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

#include <nr.h>

// Used for NR optimization: Evaluates similarity measure
extern float irtkRegistrationEvaluate(float *x);

// Used for NR optimization: Evaluates similarity measure
extern void  irtkRegistrationEvaluateDerivative(float *x, float *dx);

double irtkConjugateGradientDescentOptimizer::Run()
{
  int i, n;
  float new_similarity, old_similarity;

  // Number of variables we have to optimize
  n = _Transformation->NumberOfDOFs();

  // Current similarity
  old_similarity = _Registration->Evaluate();

  // Convert some stuff to NR
  float *x = new float[n];
  for (i = 0; i < n; i++) {
    x[i] = _Transformation->Get(i);
  }
  // Call NR Fletcher-Reeves-Polak-Ribiere routine for conjugate gradient
  frprmn(x-1, n, _Epsilon, &i, &new_similarity, irtkRegistrationEvaluate, irtkRegistrationEvaluateDerivative);

  // Convert some stuff back from NR
  for (i = 0; i < n; i++) {
    _Transformation->Put(i, x[i]);
  }

  // Delete NR memory
  delete []x;

  // Return
  if (-new_similarity > old_similarity) {
    cout << new_similarity << endl;
    return -new_similarity - old_similarity;
  } else {
    return 0;
  }
}
