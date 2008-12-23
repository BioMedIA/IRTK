/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKGENERICHISTOGRAMSIMILARITYMETRIC_H

#define _IRTKGENERICHISTOGRAMSIMILARITYMETRIC_H

/**
 * Generic class for histogram-based voxel similarity measures.
 *
 * This class implements various voxel similarity measures, including
 * sums-of-squared differences, cross-correlation, correlation ratio, mutual
 * information and normalised mutual information.
 *
 */

class irtkGenericHistogramSimilarityMetric : public irtkHistogramSimilarityMetric
{

  /// Which similarity measure are we using
  irtkSimilarityMeasure _SimilarityMeasure;

public:

  /// Constructor
  irtkGenericHistogramSimilarityMetric(irtkSimilarityMeasure = MI, int = 64, int = 64);

  /// Resets the metric.
  virtual void Reset();

  /** Resets the metric.
      \param metric To metric to reset with. */
  virtual void Reset(irtkGenericHistogramSimilarityMetric& metric);

  /** Puts the min for the metric.
      \param targetMin The target min.
      \param sourceMin The source min. */
  virtual void PutMin(double targetMin, double sourceMin);

  /** Gets the min for the metric.
      \param targetMin The target min.
      \param sourceMin The source min. */
  virtual void GetMin(double *targetMin, double *sourceMin) const;

  /** Puts the max for the metric.
      \param targetMax The target max.
      \param sourceMax The source max. */
  virtual void PutMax(double targetMax, double sourceMax);

  /** Gets the max for the metric.
      \param targetMax The target max.
      \param sourceMax The source max. */
  virtual void GetMax(double *targetMax, double *sourceMax) const;

  /** Adds to a bin for a particular time.
      \param targetBin Index to target bin.
      \param sourceBin Index to source bin.
      \param count The number to add. */
  virtual void Add(int targetBin, int sourceBin, int count = 1);

  /** Deletes from a bin for a particular time.
      \param targetBin Index to target bin.
      \param sourceBin Index to source bin.
      \param count The number to delete. */
  virtual void Delete(int targetBin, int sourceBin, int count = 1);

  /// Returns number of bins in histogram
  virtual int GetNumberOfBinsX() const;

  /// Returns number of bins in histogram
  virtual int GetNumberOfBinsY() const;

  /// Returns number of samples in histogram
  virtual int NumberOfSamples() const;

  /** Sets the similarity measure.
      \param measure The similarity measure. */
  virtual void SetMeasure(irtkSimilarityMeasure measure);

  /** Gets the similarity measure. */
  virtual irtkSimilarityMeasure GetMeasure() const;

  /// Evaluate similarity measure
  virtual double Evaluate();

};

inline irtkGenericHistogramSimilarityMetric::irtkGenericHistogramSimilarityMetric(irtkSimilarityMeasure SimilarityMeasure, int target_nbins, int source_nbins) : irtkHistogramSimilarityMetric(target_nbins, source_nbins)
{
  _SimilarityMeasure = SimilarityMeasure;
}

inline void irtkGenericHistogramSimilarityMetric::Reset()
{
  this->_histogram->Reset();
}

inline void irtkGenericHistogramSimilarityMetric::Reset(irtkGenericHistogramSimilarityMetric& metric)
{
  this->_histogram->Reset(*(metric._histogram));
}

inline void irtkGenericHistogramSimilarityMetric::PutMin(double targetMin, double sourceMin)
{
  this->_histogram->PutMin(targetMin, sourceMin);
}

inline void irtkGenericHistogramSimilarityMetric::GetMin(double *targetMin, double *sourceMin) const
{
  this->_histogram->GetMin(targetMin, sourceMin);
}

inline void irtkGenericHistogramSimilarityMetric::PutMax(double targetMax, double sourceMax)
{
  this->_histogram->PutMax(targetMax, sourceMax);
}

inline void irtkGenericHistogramSimilarityMetric::GetMax(double *targetMax, double *sourceMax) const
{
  this->_histogram->GetMax(targetMax, sourceMax);
}

inline void irtkGenericHistogramSimilarityMetric::SetMeasure(irtkSimilarityMeasure measure)
{
  _SimilarityMeasure = measure;
}

inline irtkSimilarityMeasure irtkGenericHistogramSimilarityMetric::GetMeasure() const
{
  return _SimilarityMeasure;
}

inline int irtkGenericHistogramSimilarityMetric::GetNumberOfBinsX() const
{
  return this->_histogram->GetNumberOfBinsX();
}

inline int irtkGenericHistogramSimilarityMetric::GetNumberOfBinsY() const
{
  return this->_histogram->GetNumberOfBinsY();
}

inline int irtkGenericHistogramSimilarityMetric::NumberOfSamples() const
{
  return this->_histogram->NumberOfSamples();
}

inline void irtkGenericHistogramSimilarityMetric::Add(int i, int j, int n)
{
  this->_histogram->Add(i, j, n);
}

inline void irtkGenericHistogramSimilarityMetric::Delete(int i, int j, int n)
{
  this->_histogram->Delete(i, j, n);
}

#endif
