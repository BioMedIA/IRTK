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

#define HISTORY

double irtkDownhillDescentOptimizer::Run()
{
  int i, j, k;
  double similarity, new_similarity, old_similarity;

  // Assume that the transformation is the optimal transformation
#ifdef HISTORY
  irtkHistory *history = _Registration->history;
  if (history->Find(_Transformation, similarity) == true) {
    old_similarity = new_similarity = similarity;
  } else {
    old_similarity = new_similarity = similarity = _Registration->Evaluate();
    history->Add(_Transformation, similarity);
  }
#else
  old_similarity = new_similarity = _Registration->Evaluate();
#endif
  j = 0;
  k = 0;
  // Find degree of freedom which improves the similarity measure the most
  for (i = 0; i < _Transformation->NumberOfDOFs(); i++) {
    if (_Transformation->GetStatus(i) == _Active) {
      _Transformation->Put(i, _Transformation->Get(i) + _StepSize);
#ifdef HISTORY
      if (history->Find(_Transformation, similarity) == false) {
        similarity = _Registration->Evaluate();
        history->Add(_Transformation, similarity);
      }
#else
      similarity = _Registration->Evaluate();
#endif
      if (similarity > new_similarity) {
        new_similarity = similarity;
        j =  i;
        k =  1;
      }
      _Transformation->Put(i, _Transformation->Get(i) - _StepSize);
      _Transformation->Put(i, _Transformation->Get(i) - _StepSize);
#ifdef HISTORY
      if (history->Find(_Transformation, similarity) == false) {
        similarity = _Registration->Evaluate();
        history->Add(_Transformation, similarity);
      }
#else
      similarity = _Registration->Evaluate();
#endif
      if (similarity > new_similarity) {
        new_similarity = similarity;
        j =  i;
        k = -1;
      }
      _Transformation->Put(i, _Transformation->Get(i) + _StepSize);
    }
  }
  if (new_similarity > old_similarity) {
    _Transformation->Put(j, _Transformation->Get(j) + k * _StepSize);
    return new_similarity - old_similarity;
  } else {
    return 0;
  }
}

