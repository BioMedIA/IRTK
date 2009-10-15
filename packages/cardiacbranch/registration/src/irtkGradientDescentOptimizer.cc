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

double irtkGradientDescentOptimizer::Run()
{
  int i, n;
  double similarity, new_similarity, old_similarity;

  // Assume that the transformation is the optimal transformation
  old_similarity = new_similarity = similarity = _Registration->Evaluate();

  // Number of variables we have to optimize
  n = _Transformation->NumberOfDOFs();

  // Convert some stuff to NR
  float *dx = new float[n];

  _Registration->EvaluateGradient(_StepSize, dx);

  // Step along gradient direction until no further improvement is necessary
  do {
    new_similarity = similarity;
    for (i = 0; i < _Transformation->NumberOfDOFs(); i++) {
      _Transformation->Put(i, _Transformation->Get(i) + _StepSize * dx[i]);
    }
    similarity = _Registration->Evaluate();
    if (similarity > new_similarity + _Epsilon) cout << similarity << endl;
  } while (similarity > new_similarity + _Epsilon);

  // Last step was no improvement, so back track
  for (i = 0; i < _Transformation->NumberOfDOFs(); i++) {
    _Transformation->Put(i, _Transformation->Get(i) - _StepSize * dx[i]);
  }

  // Delete NR memory
  delete []dx;

  if (new_similarity > old_similarity) {
    return new_similarity - old_similarity;
  } else {
    return 0;
  }
}
