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

#include <irtkDensity.h>

irtkDensity::irtkDensity(int nbins)
{
  int i;

  if (nbins < 1) {
    cerr << "irtkDensity::irtkDensity: Should have at least one bin";
    exit(1);
  }
  _min   = 0;
  _max   = nbins;
  _width = 1;
  _nbins = nbins;
  _norm = 0;
  if (_nbins > 0) _bins  = new double[_nbins];
  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
}

irtkDensity::irtkDensity(double min, double max, double width)
{
  int i;

  _min   = min;
  _max   = max;
  _nbins = int((_max - _min) / width);
  _width = (_max - _min) / (double)_nbins;
  _norm = 0;
  if (_nbins < 1) {
    cerr << "irtkDensity::irtkDensity: Should have at least one bin";
    exit(1);
  }
  if (_nbins > 0) _bins  = new double[_nbins];
  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
}

/*irtkDensity::irtkDensity(char *filename)
{
  int i;
  char buffer[255];

  ifstream from(filename);
  if (!from){
    cerr << "irtkDensity::Read: Can't open file " << filename << "\n";
    exit(1);
  }

  from >> buffer;
  if (strcmp(buffer, "irtkDensity") != 0){
    cerr << "irtkDensity::Read: Invalid format" << endl;
    exit(1);
  }
  from >> _nbins >> _nsamp >> _min >> _max >> _width;

  if (_nbins > 0) _bins  = new int[_nbins];
  for (i = 0; i < _nbins; i++){
    from >> _bins[i];
  }
}
*/
irtkDensity::~irtkDensity()
{
  if (_nbins > 0) {
    delete []_bins;
  }
  _nbins = 0;
  _norm = 0;
  _min   = 0;
  _max   = 0;
  _width = 0;
}

void irtkDensity::Reset()
{
  int i;

  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
  _norm = 0;
}

void irtkDensity::PutMin(double min)
{
  _min = min;
  _width = (_max - _min) / (double)_nbins;
  this->Reset();
}

void irtkDensity::PutMax(double max)
{
  _max = max;
  _width = (_max - _min) / (double)_nbins;
  this->Reset();
}

void irtkDensity::PutWidth(double width)
{
  int i;

  if (round((_max - _min) / width) < 1) {
    cerr << "irtkDensity::PutWidth: Should have at least one bin";
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
  _bins  = new double[_nbins];
  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
  _norm = 0;
}

void irtkDensity::PutNumberOfBins(int nbins)
{
  int i;

  if (nbins < 1) {
    cerr << "irtkDensity::PutNumberOfBins: Should have at least one bin";
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
  _bins  = new double[_nbins];
  for (i = 0; i < _nbins; i++) {
    _bins[i] = 0;
  }
  _norm = 0;
}

/*double irtkDensity::CDFToVal(double p)
{
  int i, sum;

  if ((p < 0) || (p > 1)){
    cerr << "irtkDensity::CDFToVal: Must be between 0 and 1" << endl;
    exit(1);
  }

  sum = 0;
  for (i = 0; i < _nbins; i++){
    sum += _bins[i];
    if (sum/(double)_nsamp >= p){
      break;
    }
  }
  return BinToVal(i);
}
*/
double irtkDensity::Mean()
{
  int i;
  double val;

  if (_norm == 0) {
    cerr << "irtkDensity::Mean: No samples in irtkDensity" << endl;
    return 0;
  }
  val = 0;
  for (i = 0; i < _nbins; i++) {
    val += _bins[i] * this->BinToVal(i);
  }
  return val / _norm;
}

double irtkDensity::Variance()
{
  int i;
  double val, mean;

  if (_norm == 0) {
    cerr << "irtkDensity::Variance: No samples in irtkDensity" << endl;
    return 0;
  }
  val  = 0;
  mean = this->Mean();
  for (i = 0; i < _nbins; i++) {
    val += _bins[i] * (this->BinToVal(i) - mean) * (this->BinToVal(i) - mean);
  }
  return val /_norm;
}

double irtkDensity::StandardDeviation()
{
  if (_norm == 0) {
    cerr << "irtkDensity::StandardDeviation: No samples in irtkDensity" << endl;
    return 0;
  }
  return sqrt(this->Variance());
}

double irtkDensity::Entropy()
{
  int i;
  double val;

  if (_norm == 0) {
    cerr << "irtkDensity::Entropy: No samples in irtkDensity" << endl;
    return 0;
  }
  val = 0;
  for (i = 0; i < _nbins; i++) {
    if (_bins[i] > 0) {
      val += (_bins[i] /_norm) * log(_bins[i] /_norm);
    }
  }
  return - val;
}

/*void irtkDensity::Read(char *filename)
{
  int i;
  char buffer[255];

  ifstream from(filename);
  if (!from){
    cerr << "irtkDensity::Read: Can't open file " << filename << "\n";
    exit(1);
  }
  if (_nbins > 0){
    delete []_bins;
    _nbins = 0;
  }

  from >> buffer;
  if (strcmp(buffer, "irtkDensity") != 0){
    cerr << "irtkDensity::Read: Invalid format" << endl;
    exit(1);
  }

  from >> _nbins >> _nsamp >> _min >> _max >> _width;
  if (_nbins > 0) _bins  = new int[_nbins];
  for (i = 0; i < _nbins; i++){
    from >> _bins[i];
  }
}

void irtkDensity::Write(char *filename)
{
  int i;

  ofstream to(filename);
  if (!to){
    cerr << "irtkDensity::Write: Can't open file " << filename << "\n";
    exit(1);
  }
  to << "irtkDensity\n";
  to << _nbins << " " << _nsamp << " " << _min << " " << _max << " "
     << _width << endl;
  for (i = 0; i < _nbins; i++){
    to << _bins[i] << endl;
  }
}
*/
void irtkDensity::Print()
{
  int i;

  cout << _nbins << " " << _norm << " " << _min << " " << _max << " "
       << _width << endl;
  for (i = 0; i < _nbins; i++) {
    cout << BinToVal(i) << ", "<< _bins[i] << endl;
  }
}

void irtkDensity::WriteToImage(char *image_name)
{

  double max_value, min_value;
  int i;
  int height = _nbins/2;

  cerr <<"Writing density to image:"<<endl;

  max_value = _bins[0];
  min_value = _bins[0];
  for (i = 1; i < _nbins; i++) {
    if (max_value < _bins[i]) max_value = _bins[i];
    if (min_value > _bins[i]) min_value = _bins[i];
  }
  cerr << "max = " << max_value <<", min = " << min_value <<endl;

  //max_value = 20849;
  //min_value = 0;

  //cerr << "max = " << max_value <<", min = " << min_value <<endl;

  cerr<<"Creating image ...";

  irtkGreyImage image(_nbins, height, 1, 1);
  image.Write(image_name);

  cerr <<"done."<<endl;

  cerr<<"Drawing density ...";
  //image.PutAsDouble( 0, 0, 0,1 );

  int old_value;

  for (i = 0; i < _nbins; i++) {

    int value = round((height-1)*(_bins[i]- min_value)/(max_value-min_value));
    if (i==0) old_value = value;
    image.Put(i, value , 0, 0, 100 );
  }

  cerr <<"done."<<endl;

  cerr<<"Writing image ..." ;
  image.Write(image_name);
  cerr <<"done."<<endl;

}
