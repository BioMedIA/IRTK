/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKLABELCONSISTENCYSIMILARITYMETRIC_H

#define _IRTKLABELCONSISTENCYSIMILARITYMETRIC_H

/**
 * Class for voxel similarity measure based on label consistency
 *
 */

class irtkLabelConsistencySimilarityMetric : public irtkHistogramSimilarityMetric
{

private:

  /// Number of samples
  double _n;

  /// Number of consistent matches
  double _match;

public:

  /// Constructor
  irtkLabelConsistencySimilarityMetric();

  /// Add sample
  virtual void Add(int, int, double =1);
  /// Remove sample
  virtual void Delete(int, int, double =1);

  /// Reset similarity metric
  virtual void Reset();

  /// Reset similarity metric
  virtual void Reset(irtkSimilarityMetric *);

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkLabelConsistencySimilarityMetric::irtkLabelConsistencySimilarityMetric()
{
  _match = 0;
  _n = 0;
}

inline void irtkLabelConsistencySimilarityMetric::Add(int x, int y, double weight)
{
  if (x == y) _match++;
  _n++;
}

inline void irtkLabelConsistencySimilarityMetric::Delete(int x, int y, double weight)
{
  if (x == y) _match--;
  _n--;
}

inline void irtkLabelConsistencySimilarityMetric::Reset()
{
  _match = 0;
  _n = 0;
}

inline void irtkLabelConsistencySimilarityMetric::Reset(irtkSimilarityMetric *metric)
{
  irtkLabelConsistencySimilarityMetric *m = dynamic_cast<irtkLabelConsistencySimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkLabelConsistencySimilarityMetric::Reset: Dynamic cast failed" << endl;
    exit(1);
  }

  _match = m->_match;
  _n = m->_n;
}

inline double irtkLabelConsistencySimilarityMetric::Evaluate()
{
  if (_n > 0) {
    return _match / _n;
  } else {
    cerr << "irtkLabelConsistencySimilarityMetric::Evaluate: No samples";
    return 0;
  }
}

#endif
