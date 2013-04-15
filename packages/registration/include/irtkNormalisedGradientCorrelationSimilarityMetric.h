/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKNORMALISEDGRADIENTCORRELATIONSIMILARITYMETRIC_H

#define _IRTKNORMALISEDGRADIENTCORRELATIONSIMILARITYMETRIC_H

/*
 * Class for voxel similarity measure based on normalised gradient correlation
 *
 * Note that the input gradient vectors must already be normalised.
 *
 */

class irtkNormalisedGradientCorrelationSimilarityMetric : public irtkSimilarityMetric
{
private:
  // Number of dimensions
  int _dim;

  // Number of samples
  double _n;

  // Measure of gradient correlation
  double _correlation;

public:
  // Constructor
  irtkNormalisedGradientCorrelationSimilarityMetric(int dim = 3);

  // Get vector dimension
  int GetDim();

  /// Add sample
  virtual void Add(int, int) {};

  /// Remove sample
  virtual void Delete(int, int) {};

  /// Add sample
  virtual void AddWeightedSample(int, int, double = 1) {};

  /// Remove sample
  virtual void DeleteWeightedSample(int, int, double = 1) {};

  // Add sample
  virtual void AddVector(double [], double [], double = 1);

  // Remove sample
  virtual void DeleteVector(double [], double [], double = 1);

  // Combine similarity metrics
  virtual void Combine(irtkSimilarityMetric *);

  // Reset similarity metric
  virtual void Reset();

  // Reset similarity metric
  virtual void ResetAndCopy(irtkSimilarityMetric *);

  // Evaluate similarity measure
  virtual double Evaluate();

};

// Get vector dimension
inline int irtkNormalisedGradientCorrelationSimilarityMetric::GetDim()
{
  return _dim;
}

inline irtkNormalisedGradientCorrelationSimilarityMetric::irtkNormalisedGradientCorrelationSimilarityMetric(int dim)
{
  _dim = dim;
  _correlation = 0;
  _n = 0;
}

inline void irtkNormalisedGradientCorrelationSimilarityMetric::AddVector(double x[], double y[], double weight)
{
  double inner_prod = 0;
  for(int i=0; i<_dim; i++){
    inner_prod += x[i] * y[i];
  }

  _correlation += inner_prod * weight;
  _n += weight;
}

inline void irtkNormalisedGradientCorrelationSimilarityMetric::DeleteVector(double x[], double y[], double weight)
{
  double inner_prod = 0;
  for(int i=0; i<_dim; i++){
    inner_prod += x[i] * y[i];
  }

  _correlation -= inner_prod * weight;
  _n -= weight;
}

inline void irtkNormalisedGradientCorrelationSimilarityMetric::Combine(irtkSimilarityMetric *metric)
{
  irtkNormalisedGradientCorrelationSimilarityMetric *m = dynamic_cast<irtkNormalisedGradientCorrelationSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkNormalisedGradientCorrelationSimilarityMetric::Combine: Dynamic cast failed" << endl;
    exit(1);
  }

  _correlation += m->_correlation;
  _n += m->_n;
}

inline void irtkNormalisedGradientCorrelationSimilarityMetric::Reset()
{
  _correlation = 0;
  _n = 0;
}

inline void irtkNormalisedGradientCorrelationSimilarityMetric::ResetAndCopy(irtkSimilarityMetric *metric)
{
  irtkNormalisedGradientCorrelationSimilarityMetric *m = dynamic_cast<irtkNormalisedGradientCorrelationSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkNormalisedGradientCorrelationSimilarityMetric::Reset: Dynamic cast failed" << endl;
    exit(1);
  }

  _dim = m->_dim;
  _correlation = m->_correlation;
  _n = m->_n;
}

inline double irtkNormalisedGradientCorrelationSimilarityMetric::Evaluate()
{
  if (_n > 0) {
    return _correlation / _n;
  } else {
    cerr << "irtkNormalisedGradientCorrelationSimilarityMetric::Evaluate: No samples";
    return 0;
  }
}

#endif
