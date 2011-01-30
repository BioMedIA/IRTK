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

#ifdef HAS_VTK

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#endif

int intersection(double x1, double y1, double x2, double y2,
                 double x3, double y3, double x4, double y4)
{
  double a = (x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3);
  double b = (x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3);
  if ((y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1) != 0) {
    a /= (y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1);
    b /= (y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1);
    if ((a >= 0) && (a < 1) && (b >= 0) && (b < 1)) {
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

void irtkPointSet::Clear()
{
  if (_m > 0) {
    delete []_data;
    _data = NULL;
    _m = 0;
    _n = 0;
  }
}

void irtkPointSet::Add(const irtkPoint &p)
{
  int i;
  irtkPoint *new_data;

  if (_n+1 <= _m) {
    // There is still enough memory left, so just add the point
    _data[_n] = p;
    _n++;
    return;
  }
  // There is not enough memory left, so allocate new point list and copy
  _m += POINTSET_SIZE;
  new_data = new irtkPoint[_m];
  for (i = 0; i < _n; i++) {
    new_data[i] = _data[i];
  }
  new_data[_n] = p;
  delete []_data;
  _data = new_data;
  _n++;
}

void irtkPointSet::Del(const irtkPoint &p)
{
  int i;
  irtkPointSet new_data(*this);

  // Delete all points
  this->Clear();
  // Copy all points except for p
  for (i = 0; i < new_data.Size(); i++) {
    if (!(new_data(i) == p)) {
      this->Add(new_data(i));
    }
  }
}

void irtkPointSet::Add(const irtkPointSet &pset)
{
  int i;

  for (i = 0; i < pset.Size(); i++) {
    this->Add(pset(i));
  }
}

void irtkPointSet::Del(const irtkPointSet &pset)
{
  int i;

  for (i = 0; i < pset.Size(); i++) {
    this->Del(pset(i));
  }
}

irtkPointSet& irtkPointSet::operator=(const irtkPointSet &pset)
{
  int i;

  // Delete old points
  this->Clear();

  // Allocate memory for new points
  if (pset._m  > 0) {
    _data = new irtkPoint[pset._m];
    _m = pset._m;
    _n = pset._n;
  }

  // Copy new points
  for (i = 0; i < _n; i++) {
    _data[i] = pset._data[i];
  }

  return *this;
}

irtkPoint irtkPointSet::CenterOfGravity() const
{
  int i;
  irtkPoint p;

  if (this->Size() == 0) {
    cerr << "irtkPointSet::CenterOfGravity(): No points in point set" << endl;
    return p;
  }
  for (i = 0; i < this->Size(); i++) {
    p += (*this)(i);
  }
  return p / (double)this->Size();
}

void irtkPointSet::BoundingBox(irtkPoint &p1, irtkPoint &p2) const
{
  int i;

  if (this->Size() == 0) {
    p1 = irtkPoint();
    p2 = irtkPoint();
    return;
  } else {
    p1 = (*this)(0);
    p2 = (*this)(0);
  }
  for (i = 1; i < this->Size(); i++) {
    p1._x = p1._x < (*this)(i)._x ? p1._x : (*this)(i)._x;
    p1._y = p1._y < (*this)(i)._y ? p1._y : (*this)(i)._y;
    p1._z = p1._z < (*this)(i)._z ? p1._z : (*this)(i)._z;
    p2._x = p2._x > (*this)(i)._x ? p2._x : (*this)(i)._x;
    p2._y = p2._y > (*this)(i)._y ? p2._y : (*this)(i)._y;
    p2._z = p2._z > (*this)(i)._z ? p2._z : (*this)(i)._z;
  }
}

ostream& operator<< (ostream& os, const irtkPointSet &pset)
{
  int i;

  os << "irtkPointSet " << pset._n << endl;
  os.setf(ios::right);
  os.setf(ios::fixed);
  os.precision(10);
  for (i = 0; i < pset._n; i++) {
    os << setw(15) << pset._data[i] << endl;
  }
  os.precision(6);
  os.unsetf(ios::right);
  os.unsetf(ios::fixed);
  return os;
}

istream& operator>> (istream& is, irtkPointSet &pset)
{
  int i, size;
  char buffer[255];
  irtkPoint p;

  // Read header
  is >> buffer;
  if (strcmp(buffer, "irtkPointSet") != 0) {
    cerr << "Can't read irtkPointSet file: " << buffer << endl;
    exit(1);
  }

  // Read size
  is >> size;

  // Read irtkPointSet
  for (i = 0; i < size; i++) {
    is >> p;
    pset.Add(p);
  }
  return is;
}

void irtkPointSet::Read(char *filename)
{
  // Clear pointset first
  this->Clear();

  // Open file stream
  ifstream from(filename);

  // Check whether file opened ok
  if (!from) {
    cerr << "irtkPointSet::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Read irtkPointSet
  from >> *this;
}

void irtkPointSet::Write(char *filename)
{
  // Open file stream
  ofstream to(filename);

  // Check whether file opened ok
  if (!to) {
    cerr << "irtkPointSet::Write: Can't open file " << filename << endl;
    exit(1);
  }

  // Write irtkPointSet
  to << *this;
}

void irtkPointSet::ReadVTK(char *filename)
{
#ifdef HAS_VTK
  int i;

  // Clear pointset first
  this->Clear();

  // Create reader
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(filename);
  reader->Update();

  // Create point set
  vtkPoints *points = reader->GetOutput()->GetPoints();

  // Convert point set
  for (i = 0; i < points->GetNumberOfPoints(); i++) {
    double p[3];
    points->GetPoint(i, p);
    this->Add(irtkPoint(p[0], p[1], p[2]));
  }
#else
  cerr << "irtkPointSet::ReadVTK: Must be compiled with VTK enabled" << endl;
  exit(1);
#endif
}

void irtkPointSet::WriteVTK(char *filename)
{
#ifdef HAS_VTK
  int i;

  // Create data
  vtkPolyData *data = vtkPolyData::New();

  // Create point set
  vtkPoints *points = vtkPoints::New();

  // Convert point set
  for (i = 0; i < this->Size(); i++) {
    double p[3];
    irtkPoint point = this->operator()(i);
    p[0] = point._x;
    p[1] = point._y;
    p[2] = point._z;
    points->InsertPoint(i, p);
  }
  data->SetPoints(points);

  // Create writer
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(filename);
  writer->SetInput(data);
  writer->Update();
#else
  cerr << "irtkPointSet::WriteVTK: Must be compiled with VTK enabled" << endl;
  exit(1);
#endif
}

irtkPoint irtkPointSet::ClosestPoint(irtkPoint& p_input){
	int i,j;
	double mindistance = 100000;
	double tmpdistance;
	irtkPoint e;

	if (this->Size() == 0) {
		cerr << "irtkPointSet::ClosestPoint(): No points in pointset" << endl;
		return e;
	}
	j = 0;
	for (i = 0; i < this->Size(); i++) {
		tmpdistance = _data[i].Distance(p_input);
		if(tmpdistance < mindistance){
			mindistance = tmpdistance;
			j = i;
		}
	}
	e = _data[j];
	return e ;
}

double irtkPointSet::PointDistance(irtkPoint& p_input){
	int i;
	double mindistance = 100000;
	double tmpdistance;

	if (this->Size() == 0) {
		cerr << "irtkPointSet::PointDistant(): No points in pointset" << endl;
		return 0;
	}
	for (i = 0; i < this->Size(); i++) {
		tmpdistance = _data[i].Distance(p_input);
		if(tmpdistance < mindistance){
			mindistance = tmpdistance;
		}
	}
	return mindistance ;
}

irtkPoint irtkPointSet::StandardDeviationEllipsoid() const
{
  int i;
  irtkPoint p;
  irtkPoint e;

  if (this->Size() == 0) {
    cerr << "irtkPointSet::StandardDeviationEllipsoid(): No points in pointset" << endl;
    return e;
  }
  p = this->CenterOfGravity();

  for (i = 0; i < this->Size(); i++) {
    e._x += ((*this)(i)._x-p._x)*((*this)(i)._x-p._x);
    e._y += ((*this)(i)._y-p._y)*((*this)(i)._y-p._y);
    e._z += ((*this)(i)._z-p._z)*((*this)(i)._z-p._z);
  }
  e._x = sqrt(e._x / ((double)this->Size()-1));
  e._y = sqrt(e._y / ((double)this->Size()-1));
  e._z = sqrt(e._z / ((double)this->Size()-1));
  return e ;
}

int irtkPointSet::IsInside(double x, double y) const
{
  int i;

  // compute no of points
  if (_n == 0) return false;

  // compute centre
  double cx = 0;
  double cy = 0;
  for (i = 0; i < _n; ++i) {
    cx += _data[i]._x;
    cy += _data[i]._y;
  }
  cx /= _n;
  cy /= _n;

  // compute a point outside the polygon.
  double ox = 0, oy = 0;
  for (i = 0; i < _n; ++i) {
    double tmp;

    tmp = _data[i]._x-cx;
    if (tmp<0) tmp = -tmp;
    if (tmp>ox) ox = tmp;

    tmp = _data[i]._y-cy;
    if (tmp<0) tmp = -tmp;
    if (tmp>oy) oy = tmp;
  }
  ox = cx + ox + oy + 1;
  oy = cy + ox + oy + 1;

  // count crossings.
  int crossings = 0;
  for (i = 0; i < _n; ++i) {
    crossings += intersection(_data[i]._x, _data[i]._y, _data[(i+1)%_n]._x,
                              _data[(i+1)%_n]._y, ox, oy, x, y);
  }

  // inside iff there was an odd number of crossings.
  return crossings % 2 != 0;
}

