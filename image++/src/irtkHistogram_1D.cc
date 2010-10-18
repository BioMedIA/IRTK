/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkImage.h>

#include <irtkHistogram.h>

template <class HistogramType> irtkHistogram_1D<HistogramType>::irtkHistogram_1D(const irtkHistogram_1D &h)
{
  int i;

  _min   = h._min;
  _max   = h._max;
  _width = h._width;
  _nbins = h._nbins;
  _nsamp = h._nsamp;
  if (_nbins > 0) _bins  = new HistogramType[_nbins];
  for (i = 0; i < _nbins; i++) {
    _bins[i] = h._bins[i];
  }
}

template <class HistogramType> irtkHistogram_1D<HistogramType>::irtkHistogram_1D(int nbins)
{
  int i;

  if (nbins < 1) {
    cerr << "irtkHistogram_1D<HistogramType>::irtkHistogram_1D: Should have at least one bin";
    exit(1);
  }
  _min   = 0;
  _max   = nbins;
  _width = 1;
  _nbins = nbins;
  _nsamp = 0;
  if (_nbins > 0) _bins  = new HistogramType[_nbins];
  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
}

template <class HistogramType> irtkHistogram_1D<HistogramType>::irtkHistogram_1D(double min, double max, double width)
{
  int i;

  _min   = min;
  _max   = max;
  _nbins = int((_max - _min) / width);
  _width = (_max - _min) / (double)_nbins;
  _nsamp = 0;
  if (_nbins < 1) {
    cerr << "irtkHistogram_1D<HistogramType>::irtkHistogram_1D: Should have at least one bin";
    exit(1);
  }
  if (_nbins > 0) _bins  = new HistogramType[_nbins];
  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
}

template <class HistogramType> irtkHistogram_1D<HistogramType>::irtkHistogram_1D(char *filename)
{
  int i;
  char buffer[255];

  ifstream from(filename);
  if (!from) {
    cerr << "irtkHistogram_1D<HistogramType>::Read: Can't open file " << filename << "\n";
    exit(1);
  }

  from >> buffer;
  if (strcmp(buffer, "irtkHistogram_1D") != 0) {
    cerr << "irtkHistogram_1D<HistogramType>::Read: Invalid format" << endl;
    exit(1);
  }
  from >> _nbins >> _nsamp >> _min >> _max >> _width;

  if (_nbins > 0) _bins  = new HistogramType[_nbins];
  for (i = 0; i < _nbins; i++) {
    from >> _bins[i];
  }
}

template <class HistogramType> irtkHistogram_1D<HistogramType>::~irtkHistogram_1D()
{
  if (_nbins > 0) {
    delete []_bins;
  }
  _nbins = 0;
  _nsamp = 0;
  _min   = 0;
  _max   = 0;
  _width = 0;
}

template <class HistogramType> void irtkHistogram_1D<HistogramType>::Reset()
{
  int i;

  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
  _nsamp = 0;
}

template <class HistogramType> void irtkHistogram_1D<HistogramType>::PutMin(double min)
{
  _min = min;
  _width = (_max - _min) / (double)_nbins;
  this->Reset();
}

template <class HistogramType> void irtkHistogram_1D<HistogramType>::PutMax(double max)
{
  _max = max;
  _width = (_max - _min) / (double)_nbins;
  this->Reset();
}

template <class HistogramType> void irtkHistogram_1D<HistogramType>::PutWidth(double width)
{
  int i;

  if (round((_max - _min) / width) < 1) {
    cerr << "irtkHistogram_1D<HistogramType>::PutWidth: Should have at least one bin";
    exit(1);
  }

  // Delete old bins if necessary
  if (_nbins > 0) {
    delete []_bins;
  }

  // Recalculate number of bins (and width)
  _nbins = round((_max - _min) / width);
  _width = (_max - _min) / (double)_nbins;

  // Allocate memory
  _bins  = new HistogramType[_nbins];
  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
  _nsamp = 0;
}

template <class HistogramType> void irtkHistogram_1D<HistogramType>::PutNumberOfBins(int nbins)
{
  int i;

  if (nbins < 1) {
    cerr << "irtkHistogram_1D<HistogramType>::PutNumberOfBins: Should have at least one bin";
    exit(1);
  }

  // Delete old bins if necessary
  if (_nbins > 0) {
    delete []_bins;
  }

  // Recalculate width of bins
  _nbins = nbins;
  _width = (_max - _min) / (double)_nbins;

  // Allocate memory
  _bins  = new HistogramType[_nbins];
  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
  _nsamp = 0;
}

template <class HistogramType> double irtkHistogram_1D<HistogramType>::CDFToVal(double p)
{
  int i;
  HistogramType sum;

  if ((p < 0) || (p > 1)) {
    cerr << "irtkHistogram_1D<HistogramType>::CDFToVal: Must be between 0 and 1" << endl;
    exit(1);
  }

  sum = 0;
  for (i = 0; i < _nbins; i++) {
    sum += _bins[i];
    if (sum/(double)_nsamp >= p) {
      break;
    }
  }
  return BinToVal(i);
}

template <> void irtkHistogram_1D<double>::Log()
{
  int i;

  if (_nsamp == 0) {
    cerr << "irtkHistogram_1D<HistogramType>::Log: No samples in irtkHistogram" << endl;
    return;
  }
  for (i = 0; i < _nbins; i++) {
    if (_bins[i] > 0) {
      _bins[i] = log(static_cast<double>(_bins[i]/(double)_nsamp));
    } else {
      _bins[i] = 0;
    }
  }
}

template <class HistogramType> double irtkHistogram_1D<HistogramType>::Mean()
{
  int i;
  double val;

  if (_nsamp == 0) {
    cerr << "irtkHistogram_1D<HistogramType>::Mean: No samples in irtkHistogram" << endl;
    return 0;
  }
  val = 0;
  for (i = 0; i < _nbins; i++) {
    val += _bins[i] * this->BinToVal(i);
  }
  return val / (double)_nsamp;
}

template <class HistogramType> double irtkHistogram_1D<HistogramType>::Variance()
{
  int i;
  double val, mean;

  if (_nsamp == 0) {
    cerr << "irtkHistogram_1D<HistogramType>::Variance: No samples in irtkHistogram" << endl;
    return 0;
  }
  val  = 0;
  mean = this->Mean();
  for (i = 0; i < _nbins; i++) {
    val += _bins[i] * (this->BinToVal(i) - mean) * (this->BinToVal(i) - mean);
  }
  return val / (double)_nsamp;
}

template <class HistogramType> double irtkHistogram_1D<HistogramType>::StandardDeviation()
{
  if (_nsamp == 0) {
    cerr << "irtkHistogram_1D<HistogramType>::StandardDeviation: No samples in irtkHistogram" << endl;
    return 0;
  }
  return sqrt(this->Variance());
}

template <class HistogramType> double irtkHistogram_1D<HistogramType>::Entropy()
{
  int i;
  double val;

  if (_nsamp == 0) {
    cerr << "irtkHistogram_1D<HistogramType>::Entropy: No samples in irtkHistogram" << endl;
    return 0;
  }
  val = 0;
  for (i = 0; i < _nbins; i++) {
    if (_bins[i] > 0) {
      val += (_bins[i] / (double)_nsamp) * log(_bins[i] / (double)_nsamp);
    }
  }
  return - val;
}

template <> void irtkHistogram_1D<double>::Smooth()
{
  int i, k;
  double *tmp, value;

  if (_nsamp == 0) {
    cerr << "irtkHistogram_1D<HistogramType>::Smooth: No samples in Histogram" << endl;
    return;
  }

  // Smoothing kernel
  double kernel[3] = { 1.0/6.0, 2.0/6.0, 1.0/6.0 };

  // Allocate temporary memory
  tmp  = new double[_nbins];

  // Smooth
  for (i = 0; i < _nbins; i++) {
    value = 0;
    for (k = 0; k < 3; k++) {
      if ((i-1+k >= 0) && (i-1+k < _nbins)) {
        value += kernel[k] * _bins[i-1+k];
      }
    }
    tmp[i] = value;
  }

  // Copy smoothed histogram back
  for (i = 0; i < _nbins; i++) _bins[i] = tmp[i];

  // Free tmp memory
  delete tmp;
}


template <class HistogramType> void irtkHistogram_1D<HistogramType>::Read(char *filename)
{
  int i;
  char buffer[255];

  ifstream from(filename);
  if (!from) {
    cerr << "irtkHistogram_1D<HistogramType>::Read: Can't open file " << filename << "\n";
    exit(1);
  }
  if (_nbins > 0) {
    delete []_bins;
    _nbins = 0;
  }

  from >> buffer;
  if (strcmp(buffer, "irtkHistogram_1D") != 0) {
    cerr << "irtkHistogram_1D<HistogramType>::Read: Invalid format" << endl;
    exit(1);
  }

  from >> _nbins >> _nsamp >> _min >> _max >> _width;
  if (_nbins > 0) _bins  = new HistogramType[_nbins];
  for (i = 0; i < _nbins; i++) {
    from >> _bins[i];
  }
}

template <class HistogramType> void irtkHistogram_1D<HistogramType>::Write(char *filename)
{
  int i;

  ofstream to(filename);
  if (!to) {
    cerr << "irtkHistogram_1D<HistogramType>::Write: Can't open file " << filename << "\n";
    exit(1);
  }
  to << "irtkHistogram_1D\n";
  to << _nbins << " " << _nsamp << " " << _min << " " << _max << " "
  << _width << endl;
  for (i = 0; i < _nbins; i++) {
    to << _bins[i] << endl;
  }
}

template <class HistogramType> void irtkHistogram_1D<HistogramType>::Print()
{
  int i;

  cout << _nbins << " " << _nsamp << " " << _min << " " << _max << " "
  << _width << endl;
  for (i = 0; i < _nbins; i++) {
    cout << _bins[i] << endl;
  }
}

template class irtkHistogram_1D<int>;
template class irtkHistogram_1D<double>;

