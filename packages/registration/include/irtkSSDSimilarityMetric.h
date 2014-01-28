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

#ifndef _IRTKSSDSIMILARITYMETRIC_H

#define _IRTKSSDSIMILARITYMETRIC_H

/**
 * Class for voxel similarity measures based on sums-of-squared differences.
 *
 */

class irtkSSDSimilarityMetric : public irtkSimilarityMetric
{

private:

  /// Number of samples
  double _n;

  /// Sums-of-squared differences
  double _ssd;

public:

  /// Constructor
  irtkSSDSimilarityMetric();

  /// Add sample
  virtual void Add(int, int);

  /// Remove sample
  virtual void Delete(int, int);

  /// Add sample
  virtual void AddWeightedSample(int, int, double = 1);

  /// Remove sample
  virtual void DeleteWeightedSample(int, int, double = 1);

  /// Combine similarity metrics
  virtual void Combine(irtkSimilarityMetric *);

  /// Reset similarity metric
  virtual void Reset();

  /// Reset similarity metric
  virtual void ResetAndCopy(irtkSimilarityMetric *);

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkSSDSimilarityMetric::irtkSSDSimilarityMetric()
{
  _ssd = 0;
  _n = 0;
}

inline void irtkSSDSimilarityMetric::Add(int x, int y)
{
  _ssd += (x-y)*(x-y);
  _n++;
}

inline void irtkSSDSimilarityMetric::Delete(int x, int y)
{
  _ssd -= (x-y)*(x-y);
  _n--;
}

inline void irtkSSDSimilarityMetric::AddWeightedSample(int x, int y, double weight)
{
  _ssd += (x-y)*(x-y)*weight;
  _n += weight;
}

inline void irtkSSDSimilarityMetric::DeleteWeightedSample(int x, int y, double weight)
{
  _ssd -= (x-y)*(x-y)*weight;
  _n -= weight;
}

inline void irtkSSDSimilarityMetric::Combine(irtkSimilarityMetric *metric)
{
  irtkSSDSimilarityMetric *m = dynamic_cast<irtkSSDSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkSSDSimilarityMetric::Combine: Dynamic cast failed" << endl;
    exit(1);
  }

  _n   += m->_n;
  _ssd += m->_ssd;
}


inline void irtkSSDSimilarityMetric::Reset()
{
  _ssd = 0;
  _n = 0;
}

inline void irtkSSDSimilarityMetric::ResetAndCopy(irtkSimilarityMetric *metric)
{
  irtkSSDSimilarityMetric *m = dynamic_cast<irtkSSDSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkSSDSimilarityMetric::ResetAndCopy: Dynamic cast failed" << endl;
    exit(1);
  }

  _ssd = m->_ssd;
  _n = m->_n;
}

inline double irtkSSDSimilarityMetric::Evaluate()
{
  if (_n > 0) {
    return - _ssd / _n;
  } else {
    cerr << "irtkSSDSimilarityMetric::Evaluate: No samples";
    return 0;
  }
}

#endif
