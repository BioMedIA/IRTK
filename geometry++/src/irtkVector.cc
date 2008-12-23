/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkGeometry.h>

irtkVector::irtkVector()
{
  _rows = 0;
  _vector = NULL;
}

irtkVector::irtkVector(int rows)
{
  _rows = 0;
  _vector = NULL;

  // Initialize vector with rows rows (sets entries to 0)
  this->Initialize(rows);
}

irtkVector::irtkVector(const irtkVector& v)
{
  int i;

  _rows = 0;
  _vector = NULL;

  // Initialize vector with v.rows rows
  this->Initialize(v._rows);

  // Copy vector
  for (i = 0; i < _rows; i++) {
    _vector[i] = v._vector[i];
  }
}

irtkVector::~irtkVector()
{
  if (_rows > 0) {
    delete []_vector;
  }
  _vector = NULL;
  _rows = 0;
}

void irtkVector::Initialize(int rows)
{
  int i;

  if (_rows > 0) {
    delete []_vector;
  }
  _vector = NULL;
  _rows = rows;
  if (_rows > 0) _vector = new double[_rows];

  // Initialize entries to 0
  for (i = 0; i < _rows; i++) {
    _vector[i] = 0.;
  }
}

ostream& operator<< (ostream& os, const irtkVector &v)
{
  // Write keyword
  os << "irtkVector " << v._rows << endl;

#ifndef WORDS_BIGENDIAN
  swap64((char *)v._vector, (char *)v._vector, v._rows);
#endif

  // Write binary data
  os.write((char *) &(v._vector[0]), v._rows*sizeof(double));

#ifndef WORDS_BIGENDIAN
  swap64((char *)v._vector, (char *)v._vector, v._rows);
#endif

  return os;
}

istream& operator>> (istream& is, irtkVector &v)
{
  int rows;
  char buffer[255];

  // Read header
  is >> buffer;
  if (strcmp(buffer, "irtkVector") != 0) {
    cerr << "irtkVector: Can't read file " << buffer << endl;
    exit(1);
  }

  // Read size
  is >> rows;

  // Allocate matrix
  v = irtkVector(rows);

  // Read header, skip comments
  is.get(buffer, 255);
  is.clear();
  is.seekg(1, ios::cur);

  // Read matrix
  is.read((char *) &(v._vector[0]), rows*sizeof(double));

#ifndef WORDS_BIGENDIAN
  swap64((char *)v._vector, (char *)v._vector, v._rows);
#endif

  return is;
}

void irtkVector::Print()
{
  int i;

  cout << "irtkVector " << _rows << endl;
  cout.setf(ios::right);
  cout.setf(ios::fixed);
  cout.precision(4);
  for (i = 0; i < _rows; i++) {
    cout << setw(15) << _vector[i] << endl;
  }
  cout.precision(6);
  cout.unsetf(ios::right);
  cout.unsetf(ios::fixed);
}

void irtkVector::Read(char *filename)
{
  // Open file stream
  ifstream from(filename, ios::in | ios::binary);

  // Check whether file opened ok
  if (!from) {
    cerr << "irtkVector::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Read vector
  from >> *this;
}

void irtkVector::Write(char *filename)
{
  // Open file stream
  ofstream to(filename, ios::out | ios::binary);

  // Check whether file opened ok
  if (!to) {
    cerr << "irtkVector::Write: Can't open file " << filename << endl;
    exit(1);
  }

  // Write vector
  to << *this;
}

#ifdef USE_VXL

#else

void irtkVector::Vector2NR(float *v) const
{
  int i;

  for (i = 0; i < _rows; i++) {
    v[i+1] = _vector[i];
  }
}

void irtkVector::NR2Vector(float *v)
{
  int i;

  for (i = 0; i < _rows; i++) {
    _vector[i] = v[i+1];
  }
}
void irtkVector::Vector2NR(double *v) const
{
  int i;

  for (i = 0; i < _rows; i++) {
    v[i+1] = _vector[i];
  }
}

void irtkVector::NR2Vector(double *v)
{
  int i;

  for (i = 0; i < _rows; i++) {
    _vector[i] = v[i+1];
  }
}

#endif
