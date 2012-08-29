/*=========================================================================

 Library   : Image Registration Toolkit (IRTK)
 Module    : $Id$
 Copyright : Imperial College, Department of Computing
 Visual Information Processing (VIP), 2008 onwards
 Date      : $Date$
 Version   : $Revision$
 Changes   : $Author$

 =========================================================================*/

#ifndef _IRTKHISTOGRAMSIMILARITYMETRIC_H

#define _IRTKHISTOGRAMSIMILARITYMETRIC_H

/**
 * Generic class for histogram-based voxel similarity measures.
 *
 * This class implements various voxel similarity measures, including
 * correlation ratio, mutual information and normalised mutual information.
 *
 */

class irtkHistogramSimilarityMetric: public irtkSimilarityMetric
{

protected:

	/// Histogram
	irtkHistogram_2D<double> *_histogram;

public:

	/// Constructor
	irtkHistogramSimilarityMetric(int = 64, int = 64);

	/// Destructor
	~irtkHistogramSimilarityMetric();

	/// Add sample
	virtual void Add(int, int);

	/// Remove sample
	virtual void Delete(int, int);

  virtual void AddWeightedSample(int, int, double = 1);

	/// Remove sample
  virtual void DeleteWeightedSample(int, int, double = 1);

	/// Combine similarity metrics
	virtual void Combine(irtkSimilarityMetric *);

	/// Reset similarity metric
	virtual void Reset();

	/// Reset similarity metric
	virtual void ResetAndCopy(irtkSimilarityMetric *);

	/// Return number of bins in X
	int NumberOfBinsX();

	/// Return number of bins in Y
	int NumberOfBinsY();

	/// Return pointer to histogram
	irtkHistogram_2D<double> *GetPointerToHistogram();

};

inline irtkHistogramSimilarityMetric::irtkHistogramSimilarityMetric(int nbins_x, int nbins_y)
{
	_histogram = new irtkHistogram_2D<double>(nbins_x, nbins_y);
}

inline irtkHistogramSimilarityMetric::~irtkHistogramSimilarityMetric()
{
	delete _histogram;
}

inline void irtkHistogramSimilarityMetric::Add(int x, int y)
{
	_histogram->Add(x, y);
}

inline void irtkHistogramSimilarityMetric::Delete(int x, int y)
{
	_histogram->Delete(x, y);
}

inline void irtkHistogramSimilarityMetric::AddWeightedSample(int x, int y, double weight)
{
  _histogram->Add(x, y, weight);
}

inline void irtkHistogramSimilarityMetric::DeleteWeightedSample(int x, int y, double weight)
{
  _histogram->Delete(x, y, weight);
}

inline void irtkHistogramSimilarityMetric::Combine(irtkSimilarityMetric *metric)
{
	int i, j;
	irtkHistogramSimilarityMetric *m = dynamic_cast<irtkHistogramSimilarityMetric *>(metric);

	if (m == NULL) {
		cerr << "irtkHistogramSimilarityMetric::Combine: Dynamic cast failed" << endl;
		exit(1);
	}

	if ((_histogram->NumberOfBinsX() != m->_histogram->NumberOfBinsX()) ||
			(_histogram->NumberOfBinsY() != m->_histogram->NumberOfBinsY())) {
		cerr << "irtkHistogramSimilarityMetric::Combine: Number of bins differs" << endl;
		exit(1);
	}
	for (j = 0; j < _histogram->NumberOfBinsY(); j++) {
		for (i = 0; i < _histogram->NumberOfBinsX(); i++) {
			_histogram->Add(i, j, m->_histogram->irtkHistogram_2D::operator()(i, j));
		}
	}
}

inline void irtkHistogramSimilarityMetric::Reset()
{
	_histogram->Reset();
}

inline void irtkHistogramSimilarityMetric::ResetAndCopy(irtkSimilarityMetric *metric)
{
	irtkHistogramSimilarityMetric *m = dynamic_cast<irtkHistogramSimilarityMetric *>(metric);

	if (m == NULL) {
		cerr << "irtkHistogramSimilarityMetric::ResetAndCopy: Dynamic cast failed" << endl;
		exit(1);
	}

	_histogram->Reset(*m->_histogram);
}

inline int irtkHistogramSimilarityMetric::NumberOfBinsX()
{
	return _histogram->NumberOfBinsX();
}

inline int irtkHistogramSimilarityMetric::NumberOfBinsY()
{
	return _histogram->NumberOfBinsY();
}

inline irtkHistogram_2D<double> * irtkHistogramSimilarityMetric::GetPointerToHistogram()
{
	return _histogram;
}

#include <irtkGenericHistogramSimilarityMetric.h>
#include <irtkMutualInformationSimilarityMetric.h>
#include <irtkNormalisedMutualInformationSimilarityMetric.h>
#include <irtkJointEntropySimilarityMetric.h>
#include <irtkCorrelationRatioXYSimilarityMetric.h>
#include <irtkCorrelationRatioYXSimilarityMetric.h>
#include <irtkKappaSimilarityMetric.h>
#include <irtkLabelConsistencySimilarityMetric.h>

#endif
