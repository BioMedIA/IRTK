/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKCROSSCORRELATIONSIMILARITYMETRIC_H

#define _IRTKCROSSCORRELATIONSIMILARITYMETRIC_H

/**
 * Class for voxel similarity measure based on cross correlation
 *
 */

class irtkCrossCorrelationSimilarityMetric : public irtkSimilarityMetric
{

private:

  /// Number of samples
  double _n;

  /// Internal variables
  double _xy, _x, _y, _x2, _y2;

public:

  /// Constructor
  irtkCrossCorrelationSimilarityMetric();

  /// Add sample
  virtual void Add(int, int, double = 1);

  /// Remove sample
  virtual void Delete(int, int, double = 1);

  /// Combine similarity metrics
  virtual void Combine(irtkSimilarityMetric *);

  /// Reset similarity metric
  virtual void Reset();

  /// Reset similarity metric
  virtual void Reset(irtkSimilarityMetric *);

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkCrossCorrelationSimilarityMetric::irtkCrossCorrelationSimilarityMetric()
{
  _xy = 0;
  _x  = 0;
  _y  = 0;
  _x2 = 0;
  _y2 = 0;
  _n  = 0;
}

inline void irtkCrossCorrelationSimilarityMetric::Add(int x, int y, double weight)
{
  _xy += weight*x*y;
  _x  += weight*x;
  _x2 += weight*x*x;
  _y  += weight*y;
  _y2 += weight*y*y;
  _n  += weight;
}

inline void irtkCrossCorrelationSimilarityMetric::Delete(int x, int y, double weight)
{
  _xy -= weight*x*y;
  _x  -= weight*x;
  _x2 -= weight*x*x;
  _y  -= weight*y;
  _y2 -= weight*y*y;
  _n  -= weight;
}

inline void irtkCrossCorrelationSimilarityMetric::Combine(irtkSimilarityMetric *metric)
{
  irtkCrossCorrelationSimilarityMetric *m = dynamic_cast<irtkCrossCorrelationSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkCrossCorrelationSimilarityMetric::Combine: Dynamic cast failed" << endl;
    exit(1);
  }

  _xy += m->_xy;
  _x2 += m->_x2;
  _y2 += m->_y2;
  _x  += m->_x;
  _y  += m->_y;
  _n  += m->_n;

}


inline void irtkCrossCorrelationSimilarityMetric::Reset()
{
  _xy = 0;
  _x2 = 0;
  _y2 = 0;
  _x  = 0;
  _y  = 0;
  _n  = 0;
}

inline void irtkCrossCorrelationSimilarityMetric::Reset(irtkSimilarityMetric *metric)
{
  irtkCrossCorrelationSimilarityMetric *m = dynamic_cast<irtkCrossCorrelationSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkCrossCorrelationSimilarityMetric::Reset: Dynamic cast failed" << endl;
    exit(1);
  }

  _xy = m->_xy;
  _x  = m->_x;
  _x2 = m->_x2;
  _y  = m->_y;
  _y2 = m->_y2;
  _n  = m->_n;
}

inline double irtkCrossCorrelationSimilarityMetric::Evaluate()
{
  if (_n > 0) {
    return (_xy - (_x * _y) / _n) / (sqrt(_x2 - _x * _x / _n) * sqrt(_y2 - _y *_y / _n));
  } else {
    cerr << "irtkCrossCorrelationSimilarityMetric::Evaluate: No samples";
    return 0;
  }
}

#endif
