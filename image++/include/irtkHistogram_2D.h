/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKHISTOGRAM_2D_H

#define _IRTKHISTOGRAM_2D_H

/** Class for 2D histograms.
 *
 *  This class defines and implements 2D histograms.
 */

template <class HistogramType> class irtkHistogram_2D : public irtkObject
{

  /// Number of bins in x-direction
  int _nbins_x;

  /// Number of bins in x-direction
  int _nbins_y;

  /// Number of samples
  HistogramType _nsamp;

  /// Min. value for x-samples, everything below is ignored
  double _min_x;

  /// Min. value for y-samples, everything below is ignored
  double _min_y;

  /// Max. value for x-samples, everything above is ignored
  double _max_x;

  /// Max. value for y-samples, everything above is ignored
  double _max_y;

  /// Width of bin in x-direction
  double _width_x;

  /// Width of bin in y-direction
  double _width_y;

  /// Dynamic memory for bins
  HistogramType **_bins;

public:

  /// Construct a histogram from another histogram
  irtkHistogram_2D(const irtkHistogram_2D &);

  /// Construct a histogram with 256 bins and samples ranging from 0 to 255
  irtkHistogram_2D(int nbins_x = 256, int nbins_y = 256);

  /// Construct a histogram for samples ranging from min to max and width
  irtkHistogram_2D(double min_x, double max_x, double width_x,
                   double min_y, double max_y, double width_y);

  /// Destructor
  ~irtkHistogram_2D(void);

  /// Clear and reset histogram
  void Reset();

  /// Clear and copy histogram
  void Reset(const irtkHistogram_2D &);

  /// Get number of bins in x-direction
  int  NumberOfBinsX() const;

  /// Put number of bins in x-direction
  void PutNumberOfBinsX(int);

  /// Get number of bins in x-direction
  int  NumberOfBinsY() const;

  /// Put number of bins in x-direction
  void PutNumberOfBinsY(int);

  /// Get number of bins in x- and y-direction
  void GetNumberOfBins(int *, int *) const;

  /// Put number of bins in x- and y-direction
  void PutNumberOfBins(int, int);

  /// Get minimum value in histogram
  void   GetMin(double *, double *) const;

  /// Put minimum value in histogram
  void   PutMin(double, double);

  /// Get maximum value in histogram
  void   GetMax(double *, double *) const;

  /// Put maximum value in histogram
  void   PutMax(double, double);

  /// Get width of bins in histogram
  void   GetWidth(double *, double *) const;

  /// Put width of bins in histogram
  void   PutWidth(double, double);

  /// Get number of samples in histogram
  HistogramType NumberOfSamples() const;

  /// Get number of samples in bin(i, j)
  HistogramType operator()(int, int) const;

  /// Add counts to bins
  void Add(int, int, HistogramType = 1);

  /// Delete counts from bins
  void Delete(int, int, HistogramType = 1);

  /// Add samples
  void AddSample(double, double, HistogramType = 1);

  /// Delete samples
  void DelSample(double, double, HistogramType = 1);

  /// Convert sample value to bin index
  int  ValToBinX(double val);

  /// Convert bin index to sample value
  double BinToValX(int    bin);

  /// Convert sample value to bin index
  int  ValToBinY(double val);

  /// Convert bin index sample value
  double BinToValY(int    bin);

  /// Convert 2D Histogram to 1D Histogram
  void HistogramX(irtkHistogram_1D<HistogramType> &);

  /// Convert 2D Histogram to 1D Histogram
  void HistogramY(irtkHistogram_1D<HistogramType> &);

  /// Log transform histogram
  void Log();

  /// Smooth histogram
  void Smooth();

  /// Calculate joint probability p(x, y)
  double JointProbability(int, int);

  /// Calculate marginal probability p(x)
  double MarginalProbabilityX(int);

  /// Calculate marginal probability p(y)
  double MarginalProbabilityY(int);

  /// Calculate conditional probability p(x|y)
  double ConditionalProbabilityXY(int, int);

  /// Calculate conditional probability p(y|x)
  double ConditionalProbabilityYX(int, int);

  /// Calculate mean
  double MeanX();

  /// Calculate mean
  double MeanY();

  /// Calculate conditional mean
  double ConditionalMeanXY(int);

  /// Calculate conditional mean
  double ConditionalMeanYX(int);

  /// Calculate variance
  double VarianceX();

  /// Calculate variance
  double VarianceY();

  /// Calculate conditional variance
  double ConditionalVarianceXY(int);

  /// Calculate conditional variance
  double ConditionalVarianceYX(int);

  /// Calculate covariance
  double Covariance();

  /// Calculate standard deviation
  double StandardDeviationX();

  /// Calculate standard deviation
  double StandardDeviationY();

  /// Calculate marginal entropy
  double EntropyX();

  /// Calculate marginal entropy
  double EntropyY();

  /// Calculate joint entropy
  double JointEntropy();

  /// Calculate mutual information
  double MutualInformation();

  /// Calculate normalized mutual information
  double NormalizedMutualInformation();

  /// Calculate cross correlation
  double CrossCorrelation();

  /// Calculate correlation ratio
  double CorrelationRatioXY();

  /// Calculate correlation ratio
  double CorrelationRatioYX();

  /// Calculate sums of squared differences
  double SumsOfSquaredDifferences();

  /// Calcualate label consistency
  double LabelConsistency();

  /// Calcualate kappa statistic
  double Kappa();

  /// Read histogram
  void Read (char *filename);

  /// Write histogram
  void Write(char *filename);

  /// Write histogram as 2D image
  void WriteAsImage(char *filename);

  /// Print histogram
  void Print();
};

template <class HistogramType> inline int irtkHistogram_2D<HistogramType>::NumberOfBinsX() const
{
  return _nbins_x;
}

template <class HistogramType> inline int irtkHistogram_2D<HistogramType>::NumberOfBinsY() const
{
  return _nbins_y;
}

template <class HistogramType> inline HistogramType irtkHistogram_2D<HistogramType>::NumberOfSamples() const
{
  return _nsamp;
}

template <class HistogramType> inline HistogramType irtkHistogram_2D<HistogramType>::operator()(int i, int j) const
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins_x) || (j < 0) || (j >= _nbins_y)) {
    cerr << "irtkHistogram_2D<HistogramType>::operator(): No such bin" << endl;
    exit(1);
  }
#endif
  return _bins[j][i];
}

template <class HistogramType> inline void irtkHistogram_2D<HistogramType>::Add(int i, int j, HistogramType n)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins_x) || (j < 0) || (j >= _nbins_y)) {
    cerr << "irtkHistogram_2D<HistogramType>::Add: No such bin " << i << " " << j << endl;
    exit(1);
  }
#endif
  _bins[j][i] += n;
  _nsamp      += n;
}

template <class HistogramType> inline void irtkHistogram_2D<HistogramType>::Delete(int i, int j, HistogramType n)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins_x) || (j < 0) || (j >= _nbins_y)) {
    cerr << "irtkHistogram_2D<HistogramType>::Delete: No such bin " << i << " " << j << endl;
    exit(1);
  }
#endif
  _bins[j][i] -= n;
  _nsamp      -= n;
}

template <class HistogramType> inline int irtkHistogram_2D<HistogramType>::ValToBinX(double val)
{
  int index;

#ifndef NO_BOUNDS
  if ((val < _min_x) || (val > _max_x)) {
    cerr << "irtkHistogram_2D<HistogramType>::ValToBinX: Must be between " << _min_x << " and "
         << _max_x << endl;
    exit(1);
  }
#endif
  index = round(_nbins_x * (val - _min_x - 0.5*_width_x) / (_max_x - _min_x));
  if (index < 0) index = 0;
  if (index > _nbins_x-1) index = _nbins_x - 1;
  return index;
}

template <class HistogramType> inline int irtkHistogram_2D<HistogramType>::ValToBinY(double val)
{
  int index;

#ifndef NO_BOUNDS
  if ((val < _min_y) || (val > _max_y)) {
    cerr << "irtkHistogram_2D<HistogramType>::ValToBinY: Must be between " << _min_y << " and "
         << _max_y << endl;
    exit(1);
  }
#endif
  index = round(_nbins_y * (val - _min_y - 0.5*_width_y) / (_max_y - _min_y));
  if (index < 0) index = 0;
  if (index > _nbins_y-1) index = _nbins_y - 1;
  return index;
}

template <class HistogramType> inline double irtkHistogram_2D<HistogramType>::BinToValX(int i)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins_x)) {
    cerr << "irtkHistogram_1D::BinToValX: Must be between 0 and " << _nbins_x << endl;
    exit(1);
  }
#endif
  return (i*(_max_x - _min_x)/_nbins_x + _min_x) + _width_x/2.0;
}

template <class HistogramType> inline double irtkHistogram_2D<HistogramType>::BinToValY(int i)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins_y)) {
    cerr << "irtkHistogram_1D::BinToValY: Must be between 0 and " << _nbins_y << endl;
    exit(1);
  }
#endif
  return (i*(_max_y - _min_y)/_nbins_y + _min_y) + _width_y/2.0;
}

template <class HistogramType> inline double irtkHistogram_2D<HistogramType>::JointProbability(int i, int j)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins_x) || (j < 0) || (j >= _nbins_y)) {
    cerr << "irtkHistogram_1D::JointProbability: No such bin " << i << " ";
    cerr << j << endl;
    exit(1);
  }
#endif
  return _bins[j][i] / (double) _nsamp;
}

template <class HistogramType> inline double irtkHistogram_2D<HistogramType>::MarginalProbabilityX(int i)
{
  int j, n;

#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins_x)) {
    cerr << "irtkHistogram_1D::MarginalProbabilityX: No such bin " << i << endl;
    exit(1);
  }
#endif
  n = 0;
  for (j = 0; j < _nbins_y; j++) {
    n += _bins[j][i];
  }
  return n / (double) _nsamp;
}

template <class HistogramType> inline double irtkHistogram_2D<HistogramType>::MarginalProbabilityY(int i)
{
  int j, n;

#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins_y)) {
    cerr << "irtkHistogram_1D::MarginalProbabilityY: No such bin " << i << endl;
    exit(1);
  }
#endif
  n = 0;
  for (j = 0; j < _nbins_x; j++) {
    n += _bins[i][j];
  }
  return n / (double) _nsamp;
}

template <class HistogramType> inline double irtkHistogram_2D<HistogramType>::ConditionalProbabilityXY(int i, int j)
{
  double p;

  p = this->MarginalProbabilityY(j);
  if (p > 0) {
    return this->JointProbability(i, j) / p;
  } else {
    return 0;
  }
}

template <class HistogramType> inline double irtkHistogram_2D<HistogramType>::ConditionalProbabilityYX(int i, int j)
{
  double p;

  p = this->MarginalProbabilityX(j);
  if (p > 0) {
    return this->JointProbability(j, i) / p;
  } else {
    return 0;
  }
}

#endif
