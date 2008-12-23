/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKDENSITY_H

#define _IRTKDENSITY_H

/** Class for non-parametric density.
 *
 *  This class defines and implements non-parametric density.
 */

class irtkDensity : public irtkObject
{

  /// Number of bins
  int _nbins;

  /// Sum of all values (integral)
  double _norm;

  /// Min. value for samples, everything below is ignored
  double _min;

  /// Max. value for samples, everything above is ignored
  double _max;

  /// Width of bins
  double _width;

  /// Dynamic memory for bins
  double *_bins;

public:

  /// Construct a histogram from another histogram
//  irtkHistogram_1D(const irtkHistogram_1D &);

  /// Construct a density with 256 bins and samples ranging from 0 to 255
  irtkDensity(int nbins = 256);

  /// Construct a density for samples ranging from min to max and width
  irtkDensity(double min, double max, double width);

  /// Read constructor
  //irtkHistogram_1D(char *);

  /// Destructor
  ~irtkDensity(void);

  /// Clear and reset
  void Reset();

  /// Get number of bins
  int  GetNumberOfBins() const;

  /// Put number of bins
  void PutNumberOfBins(int);

  /// Get minimum value
  double GetMin() const;

  /// Put minimum value
  void   PutMin(double);

  /// Get maximum value
  double GetMax() const;

  /// Put maximum value
  void   PutMax(double);

  /// Get width of bins
  double GetWidth() const;

  /// Put width of bins
  void   PutWidth(double);

  /// Get norm
  double Norm() const;

  /// Get value in bin(i)
  double operator()(int) const;

  /// Add value to bin
  void Add(int, double = 1);

  /// Delete value from bin
  void Delete(int, double = 1);

  /// Add sample to bin
  void AddSample(double, double = 1);

  /// Delete sample from bin
  void DelSample(double, double = 1);

  /// Convert sample value to bin index
  int    ValToBin(double val);

  /// Convert bin index to sample value
  double BinToVal(int    bin);

  /// Convert bin into probability density distributions
  double BinToPDF(int bin);

  /// Convert sample value into probability density distributions
  double ValToPDF(double val);

  /// Convert bin into cumulative density distributions
//  double BinToCDF(int bin);

  /// Convert sample value into cumulative  density distributions
//  double ValToCDF(double val);

  /// Convert cumulative density distributions to bin value
//  double CDFToBin(double p);

  /// Convert cumulative density distributions to sample value
//  double CDFToVal(double p);

  /// Calculate mean
  double Mean();

  /// Calculate variance
  double Variance();

  /// Calculate standard deviation
  double StandardDeviation();

  /// Calculate entropy
  double Entropy();

  /// Read histogram
//  void Read (char *filename);

  /// Wrirte histogram
//  void Write(char *filename);

  /// Print histogram
  void Print();

  void WriteToImage(char *image_name);
};

inline int irtkDensity::GetNumberOfBins() const
{
  return _nbins;
}

inline double irtkDensity::GetMin() const
{
  return _min;
}

inline double irtkDensity::GetMax() const
{
  return _max;
}

inline double irtkDensity::GetWidth() const
{
  return _width;
}

inline double irtkDensity::Norm() const
{
  return _norm;
}

inline double irtkDensity::operator()(int i) const
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkDensity::operator(): No such bin" << endl;
    exit(1);
  }
#endif
  return _bins[i];
}

inline void irtkDensity::Add(int i, double value)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkDensity::Add: No such bin" << endl;
    exit(1);
  }
#endif
  _bins[i] += value;
  _norm   += value;
}

inline void irtkDensity::Delete(int i, double value)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkDensity::Delete: No such bin" << endl;
    exit(1);
  }
#endif
  _bins[i] -= value;
  _norm   -= value;

}

inline void irtkDensity::AddSample(double x, double value)
{
  if (x < _min) return;
  if (x > _max) return;
  double bin = (x - _min) / (_max - _min) *_nbins -0.5;
  /*double bin = (x - _min) / (_max - _min) *(_nbins-1) ;*/
  int index1 = (int)floor(bin);
  int index2 = (int)ceil(bin);
  if (index2 == 0) _bins[0]+=value;
  if (index1 == _nbins-1) _bins[_nbins-1]+=value;
  if ((index1>=0)&&(index2<=_nbins-1)) {
    if (index1 ==index2) _bins[index1]+=value;
    else {
      _bins[index1]+=(index2-bin)*value;
      _bins[index2]+=(bin-index1)*value;
    }
  }
  _norm += value;
}

inline void irtkDensity::DelSample(double x, double value)
{
  int index;

  if (x < _min) return;
  if (x > _max) return;
  index = ValToBin(x);
  _bins[index] -= value;
  _norm -= value;
}

inline int irtkDensity::ValToBin(double val)
{
  int index;

#ifndef NO_BOUNDS
  if ((val < _min) || (val > _max)) {
    cerr << "irtkDensity::ValToBin: Must be between " << _min << " and "
         << _max << endl;
    exit(1);
  }
#endif
  index = round((val - _min) / (_max - _min) * _nbins -0.5);
  if (index < 0) index = 0;
  if (index > _nbins-1) index = _nbins - 1;
  return index;
}

inline double irtkDensity::BinToVal(int i)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkDensity::BinToVal: Must be between 0 and " << _nbins << endl;
    exit(1);
  }
#endif
  return (_min + (i + 0.5)*(_max - _min)/_nbins);
}

inline double irtkDensity::BinToPDF(int i)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkDensity::BinToPDF: No such bin" << endl;
    exit(1);
  }
#endif
  if (_norm == 0) return 0;
  return _bins[i] / _norm;
}

inline double irtkDensity::ValToPDF(double value)
{
  int i = ValToBin(value);
  if (_norm == 0) return 0;
  return _bins[i] / _norm;
}

#endif

