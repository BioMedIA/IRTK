/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKMLSIMILARITYMETRIC_H

#define _IRTKMLSIMILARITYMETRIC_H

#include <irtkEMClassificationTemplateBiasCorrection.h>

/**
 * Class for voxel similarity measures based on maximum likelihood.
 *
 */

class irtkMLSimilarityMetric : public irtkSimilarityMetric
{

private:

  /// log likelihood
  double _ll;
  ///number of samples
  int _n;

public:
  ///
  irtkEMClassification *_classification;


public:

  /// Constructor
  irtkMLSimilarityMetric(irtkEMClassification *);

  /// Add sample
  virtual void Add(int, int);

  /// Remove sample
  virtual void Delete(int, int);

  /// Combine similarity metrics
  virtual void Combine(irtkSimilarityMetric *);

  /// Reset similarity metric
  virtual void Reset();

  /// Reset similarity metric
  virtual void ResetAndCopy(irtkSimilarityMetric *);

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkMLSimilarityMetric::irtkMLSimilarityMetric(irtkEMClassification *classification)
{
  _ll = 0;
  _n=0;
  cerr<<"member";
  _classification=classification;
  cerr<<"GInit";
  _classification->GInit();
  cerr<<"done."<<endl;

}


inline void irtkMLSimilarityMetric::Reset()
{
  _ll = 0;
  _n=0;
}

inline void irtkMLSimilarityMetric::Add(int x, int y)
{
  //cerr<<"Add"<<endl;
  //if ((x>0)&&(y>0))
  _ll += _classification->PointLogLikelihoodGMM(x,y);
  //_ll += (y-x)*(y-x);
  _n ++;
}

inline void irtkMLSimilarityMetric::Delete(int x, int y)
{
  _ll -= _classification->PointLogLikelihoodGMM(x,y);
  //  _ll -= (y-x)*(y-x);
  _n--;
}

inline void irtkMLSimilarityMetric::Combine(irtkSimilarityMetric *)
{
  cerr << "irtkMLSimilarityMetric::Combine: Not implemented" << endl;
  exit(1);
}

inline void irtkMLSimilarityMetric::ResetAndCopy(irtkSimilarityMetric *metric)
{
  irtkMLSimilarityMetric *m = dynamic_cast<irtkMLSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkMLSimilarityMetric::ResetAndCopy: Dynamic cast failed" << endl;
    exit(1);
  }

  _ll = m->_ll;
  _classification = m->_classification;
  _n=m->_n;
}

inline double irtkMLSimilarityMetric::Evaluate()
{
  return -_ll/_n;
}

#endif
