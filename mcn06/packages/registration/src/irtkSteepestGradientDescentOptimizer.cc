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

double irtkSteepestGradientDescentOptimizer::Run()
{
  int i, n;
  float new_similarity, old_similarity;

  // Number of variables we have to optimize
  n = _Transformation->NumberOfDOFs();

  // Current similarity
  old_similarity = _Registration->Evaluate();
  cout << old_similarity << endl;

  // Convert some stuff to NR
  float *x  = new float[n];
  float *dx = new float[n];
  for (i = 0; i < n; i++) {
    x[i] = _Transformation->Get(i);
  }

  _Registration->EvaluateGradient(_StepSize, dx);

  // Do line minimization
  linmin(x-1, dx-1, n, &new_similarity, irtkRegistrationEvaluate);

  // Convert some stuff back from NR
  for (i = 0; i < n; i++) {
    _Transformation->Put(i, x[i]);
  }
  cout << -new_similarity << " " << endl;

  // Delete NR memory
  delete []x;
  delete []dx;

  // Return
  if (-new_similarity > old_similarity) {
    return -new_similarity - old_similarity;
  } else {
    return 0;
  }
}
