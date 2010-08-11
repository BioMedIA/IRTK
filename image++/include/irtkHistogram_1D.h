/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKHISTOGRAM_1D_H

#define _IRTKHISTOGRAM_1D_H

/** Class for 2D histograms.
 *
 *  This class defines and implements 2D histograms.
 */

template <class HistogramType> class irtkHistogram_1D : public irtkObject
{

  /// Number of bins
  int _nbins;

  /// Number of samples
  HistogramType _nsamp;

  /// Min. value for samples, everything below is ignored
  double _min;

  /// Max. value for samples, everything below is ignored
  double _max;

  /// Width of bins
  double _width;

  /// Dynamic memory for bins
  HistogramType *_bins;

public:

  /// Construct a histogram from another histogram
  irtkHistogram_1D(const irtkHistogram_1D &);

  /// Construct a histogram with 256 bins and samples ranging from 0 to 255
  irtkHistogram_1D(int nbins = 256);

  /// Construct a histogram for samples ranging from min to max and width
  irtkHistogram_1D(double min, double max, double width);

  /// Read constructor
  irtkHistogram_1D(char *);

  /// Destructor
  ~irtkHistogram_1D(void);

  /// Clear and reset histogram
  void Reset();

  /// Get number of bins in histogram
  int  NumberOfBins() const;

  /// Put number of bins in histogram
  void PutNumberOfBins(int);

  /// Get minimum value in histogram
  double GetMin() const;

  /// Put minimum value in histogram
  void   PutMin(double);

  /// Get maximum value in histogram
  double GetMax() const;

  /// Put maximum value in histogram
  void   PutMax(double);

  /// Get width of bins in histogram
  double GetWidth() const;

  /// Put width of bins in histogram
  void   PutWidth(double);

  /// Get number of samples in histogram
  HistogramType NumberOfSamples() const;

  /// Get number of samples in bin(i)
  HistogramType operator()(int) const;

  /// Add counts to bin
  void Add(int, HistogramType = 1);

  /// Delete counts from bin
  void Delete(int, HistogramType = 1);

  /// Add sample to bin
  void AddSample(double, HistogramType = 1);

  /// Delete sample from bin
  void DelSample(double, HistogramType = 1);

  /// Convert sample value to bin index
  int    ValToBin(double val);

  /// Convert bin index to sample value
  double BinToVal(int    bin);

  /// Convert bin into probability density distributions
  double BinToPDF(int bin);

  /// Convert sample value into probability density distributions
  double ValToPDF(double val);

  /// Convert bin into cumulative density distributions
  double BinToCDF(int bin);

  /// Convert sample value into cumulative  density distributions
  double ValToCDF(double val);

  /// Convert cumulative density distributions to bin value
  double CDFToBin(double p);

  /// Convert cumulative density distributions to sample value
  double CDFToVal(double p);

  /// Log transform histogram
  void Log();

  /// Smooth histogram
  void Smooth();

  /// Calculate mean
  double Mean();

  /// Calculate variance
  double Variance();

  /// Calculate standard deviation
  double StandardDeviation();

  /// Calculate entropy
  double Entropy();

  /// Read histogram
  void Read (char *filename);

  /// Wrirte histogram
  void Write(char *filename);

  /// Print histogram
  void Print();
};

template <class HistogramType> inline int irtkHistogram_1D<HistogramType>::NumberOfBins() const
{
  return _nbins;
}

template <class HistogramType> inline double irtkHistogram_1D<HistogramType>::GetMin() const
{
  return _min;
}

template <class HistogramType> inline double irtkHistogram_1D<HistogramType>::GetMax() const
{
  return _max;
}

template <class HistogramType> inline double irtkHistogram_1D<HistogramType>::GetWidth() const
{
  return _width;
}

template <class HistogramType> inline HistogramType irtkHistogram_1D<HistogramType>::NumberOfSamples() const
{
  return _nsamp;
}

template <class HistogramType> inline HistogramType irtkHistogram_1D<HistogramType>::operator()(int i) const
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::operator(): No such bin" << endl;
    exit(1);
  }
#endif
  return _bins[i];
}

template <class HistogramType> inline void irtkHistogram_1D<HistogramType>::Add(int i, HistogramType n)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::Add: No such bin" << endl;
    exit(1);
  }
#endif
  _bins[i] += n;
  _nsamp   += n;
}

template <class HistogramType> inline void irtkHistogram_1D<HistogramType>::Delete(int i, HistogramType n)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::Delete: No such bin" << endl;
    exit(1);
  }
#endif
  _bins[i] -= n;
  _nsamp   -= n;
}

template <class HistogramType> inline void irtkHistogram_1D<HistogramType>::AddSample(double x, HistogramType n)
{
  int index;

  if (x < _min) return;
  if (x > _max) return;
  index = round(_nbins * (x - _min - 0.5*_width) / (_max - _min));
  if (index < 0) index = 0;
  if (index > _nbins-1) index = _nbins - 1;
  _bins[index] += n;
  _nsamp       += n;
}

template <class HistogramType> inline void irtkHistogram_1D<HistogramType>::DelSample(double x, HistogramType n)
{
  int index;

  if (x < _min) return;
  if (x > _max) return;
  index = round(_nbins * (x - _min - 0.5*_width) / (_max - _min));
  if (index < 0) index = 0;
  if (index > _nbins-1) index = _nbins - 1;
  _bins[index] -= n;
  _nsamp       -= n;
}

template <class HistogramType> inline int irtkHistogram_1D<HistogramType>::ValToBin(double val)
{
  int index;

#ifndef NO_BOUNDS
  if ((val < _min) || (val > _max)) {
    cerr << "irtkHistogram_1D<HistogramType>::ValToBin: Must be between " << _min << " and "
         << _max << endl;
    exit(1);
  }
#endif
  index = round(_nbins * (val - _min - 0.5*_width) / (_max - _min));
  if (index < 0) index = 0;
  if (index > _nbins-1) index = _nbins - 1;
  return index;
}

template <class HistogramType> inline double irtkHistogram_1D<HistogramType>::BinToVal(int i)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i > _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::BinToVal: Must be between 0 and " << _nbins << endl;
    exit(1);
  }
#endif
  return (i*(_max - _min)/_nbins + _min) + _width/2.0;
}

template <class HistogramType> inline double irtkHistogram_1D<HistogramType>::BinToPDF(int i)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::BinToPDF: No such bin" << endl;
    exit(1);
  }
#endif
  if (_nsamp == 0) return 0;
  return _bins[i] / (double)_nsamp;
}

#endif
