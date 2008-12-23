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

class irtkHistogram_1D : public irtkObject
{

  /// Number of bins
  int _nbins;

  /// Number of samples
  int _nsamp;

  /// Min. value for samples, everything below is ignored
  double _min;

  /// Max. value for samples, everything below is ignored
  double _max;

  /// Width of bins
  double _width;

  /// Dynamic memory for bins
  int *_bins;

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
  int  GetNumberOfBins() const;

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
  int NumberOfSamples() const;

  /// Get number of samples in bin(i)
  int operator()(int) const;

  /// Add counts to bin
  void Add(int, int = 1);

  /// Delete counts from bin
  void Delete(int, int = 1);

  /// Add sample to bin
  void AddSample(double);

  /// Delete sample from bin
  void DelSample(double);

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

inline int irtkHistogram_1D::GetNumberOfBins() const
{
  return _nbins;
}

inline double irtkHistogram_1D::GetMin() const
{
  return _min;
}

inline double irtkHistogram_1D::GetMax() const
{
  return _max;
}

inline double irtkHistogram_1D::GetWidth() const
{
  return _width;
}

inline int irtkHistogram_1D::NumberOfSamples() const
{
  return _nsamp;
}

inline int irtkHistogram_1D::operator()(int i) const
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D::operator(): No such bin" << endl;
    exit(1);
  }
#endif
  return _bins[i];
}

inline void irtkHistogram_1D::Add(int i, int n)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D::Add: No such bin" << endl;
    exit(1);
  }
#endif
  _bins[i] += n;
  _nsamp   += n;
}

inline void irtkHistogram_1D::Delete(int i, int n)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D::Delete: No such bin" << endl;
    exit(1);
  }
#endif
  _bins[i] -= n;
  _nsamp   -= n;
}

inline void irtkHistogram_1D::AddSample(double x)
{
  int index;

  if (x < _min) return;
  if (x > _max) return;
  index = round(_nbins * (x - _min - 0.5*_width) / (_max - _min));
  if (index < 0) index = 0;
  if (index > _nbins-1) index = _nbins - 1;
  _bins[index]++;
  _nsamp++;
}

inline void irtkHistogram_1D::DelSample(double x)
{
  int index;

  if (x < _min) return;
  if (x > _max) return;
  index = round(_nbins * (x - _min - 0.5*_width) / (_max - _min));
  if (index < 0) index = 0;
  if (index > _nbins-1) index = _nbins - 1;
  _bins[index]--;
  _nsamp--;
}

inline int irtkHistogram_1D::ValToBin(double val)
{
  int index;

#ifndef NO_BOUNDS
  if ((val < _min) || (val > _max)) {
    cerr << "irtkHistogram_1D::ValToBin: Must be between " << _min << " and "
         << _max << endl;
    exit(1);
  }
#endif
  index = round(_nbins * (val - _min - 0.5*_width) / (_max - _min));
  if (index < 0) index = 0;
  if (index > _nbins-1) index = _nbins - 1;
  return index;
}

inline double irtkHistogram_1D::BinToVal(int i)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i > _nbins)) {
    cerr << "irtkHistogram_1D::BinToVal: Must be between 0 and " << _nbins << endl;
    exit(1);
  }
#endif
  return (i*(_max - _min)/_nbins + _min) + _width/2.0;
}

inline double irtkHistogram_1D::BinToPDF(int i)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D::BinToPDF: No such bin" << endl;
    exit(1);
  }
#endif
  if (_nsamp == 0) return 0;
  return _bins[i] / (double)_nsamp;
}

#endif
