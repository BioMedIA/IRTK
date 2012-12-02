/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkTransformation.h>

#include <newt2.h>

void irtkFreeFormTransformation4D::UpdateMatrix()
{
  // Update image to world coordinate system matrix
  _matL2W.Ident();
  _matL2W(0, 0) = _xaxis[0];
  _matL2W(1, 0) = _xaxis[1];
  _matL2W(2, 0) = _xaxis[2];
  _matL2W(0, 1) = _yaxis[0];
  _matL2W(1, 1) = _yaxis[1];
  _matL2W(2, 1) = _yaxis[2];
  _matL2W(0, 2) = _zaxis[0];
  _matL2W(1, 2) = _zaxis[1];
  _matL2W(2, 2) = _zaxis[2];

  irtkMatrix tmp1(4, 4);
  tmp1.Ident();
  tmp1(0, 3) = - (_x - 1) / 2.0;
  tmp1(1, 3) = - (_y - 1) / 2.0;
  tmp1(2, 3) = - (_z - 1) / 2.0;

  irtkMatrix tmp2(4, 4);
  tmp2.Ident();
  tmp2(0, 0) = _dx;
  tmp2(1, 1) = _dy;
  tmp2(2, 2) = _dz;

  irtkMatrix tmp3(4, 4);
  tmp3.Ident();
  tmp3(0, 3) = _origin._x;
  tmp3(1, 3) = _origin._y;
  tmp3(2, 3) = _origin._z;

  _matL2W = tmp3 * (_matL2W * (tmp2 * tmp1));

  // Update world to image coordinate system matrix
  _matW2L.Ident();
  _matW2L(0, 0) = _xaxis[0];
  _matW2L(0, 1) = _xaxis[1];
  _matW2L(0, 2) = _xaxis[2];
  _matW2L(1, 0) = _yaxis[0];
  _matW2L(1, 1) = _yaxis[1];
  _matW2L(1, 2) = _yaxis[2];
  _matW2L(2, 0) = _zaxis[0];
  _matW2L(2, 1) = _zaxis[1];
  _matW2L(2, 2) = _zaxis[2];

  tmp1.Ident();
  tmp1(0, 3) = (_x - 1) / 2.0;
  tmp1(1, 3) = (_y - 1) / 2.0;
  tmp1(2, 3) = (_z - 1) / 2.0;

  tmp2.Ident();
  tmp2(0, 0) = 1.0 / _dx;
  tmp2(1, 1) = 1.0 / _dy;
  tmp2(2, 2) = 1.0 / _dz;

  tmp3.Ident();
  tmp3(0, 3) = - _origin._x;
  tmp3(1, 3) = - _origin._y;
  tmp3(2, 3) = - _origin._z;

  _matW2L = tmp1 * (tmp2 * (_matW2L * tmp3));
}

double ****irtkFreeFormTransformation4D::Allocate(double ****data, int x, int y, int z, int t)
{
  int i, j, k, l;

  if ((x == 0) || (y == 0) || (z == 0) || (t == 0)) {
    return NULL;
  }

  if ((data = new double ***[t+8]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x " << z << " x " << t << "\n";
    exit(1);
  }
  data += 4;

  if ((data[-4] = new double **[(t+8)*(z+8)]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x " << z << " x " << t << "\n";
    exit(1);
  }
  data[-4] += 4;

  for (i = -3; i < t+4; i++) {
    data[i] = data[i-1] + (z+8);
  }

  if ((data[-4][-4] = new double *[(t+8)*(z+8)*(y+8)]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x " << z << " x " << t << "\n";
    exit(1);
  }
  data[-4][-4] += 4;

  for (i = -4; i < t+4; i++) {
    for (j = -4; j < z+4; j++) {
      data[i][j] = data[-4][-4] + (i+4)*(z+8)*(y+8) + (j+4)*(y+8);
    }
  }

  if ((data[-4][-4][-4] = new double [(t+8)*(z+8)*(y+8)*(x+8)]) == NULL) {
    cerr << "Allocate: malloc failed for " << x << " x " << y << " x " << z << " x " << t << "\n";
    exit(1);
  }
  data[-4][-4][-4] += 4;

  for (i = -4; i < t+4; i++) {
    for (j = -4; j < z+4; j++) {
      for (k = -4; k < y+4; k++) {
        data[i][j][k] = data[-4][-4][-4] + (i+4)*(z+8)*(y+8)*(x+8) + (j+4)*(y+8)*(x+8) + (k+4)*(x+8);
      }
    }
  }

  for (i = -4; i < t+4; i++) {
    for (j = -4; j < z+4; j++) {
      for (k = -4; k < y+4; k++) {
        for (l = -4; l < x+4; l++) {
          data[i][j][k][l] = 0;
        }
      }
    }
  }
  return data;
}

double ****irtkFreeFormTransformation4D::Deallocate(double ****data, int, int, int, int)
{
  if (data != NULL) {
    delete [](data[-4][-4][-4]-4);
    delete [](data[-4][-4]-4);
    delete [](data[-4]-4);
    delete [](data-4);
  }
  return NULL;
}

bool irtkFreeFormTransformation4D::IsIdentity()
{
  int i, j, k;

  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        if (_xdata[k][j][i] != 0) return false;
        if (_ydata[k][j][i] != 0) return false;
        if (_zdata[k][j][i] != 0) return false;
      }
    }
  }
  return true;
}

void irtkFreeFormTransformation4D::PutStatusCP(int i, int j, int k, int l, _Status status_x, _Status status_y, _Status status_z)
{
  if ((i < 0) || (i >= _x) || (j < 0) || (j >= _y) || (k < 0) || (k >= _z) || (l < 0) || (l >= _t)) {
    cerr << "irtkFreeFormTransformation4D::PutStatus: No such control point"
         << endl;
    exit(1);
  }

  int index = this->LatticeToIndex(i, j, k, l);
  _status[index] = status_x;
  _status[index+_x*_y*_z] = status_y;
  _status[index+2*_x*_y*_z] = status_z;
}

void irtkFreeFormTransformation4D::GetStatusCP(int i, int j, int k, int l,  _Status &status_x, _Status &status_y, _Status &status_z)
{
  if ((i < 0) || (i >= _x) || (j < 0) || (j >= _y) || (k < 0) || (k >= _z) || (l < 0) || (l >= _t)) {
    cerr << "irtkFreeFormTransformation4D::GetStatus: No such control point"
         << endl;
    exit(1);
  }

  int index = this->LatticeToIndex(i, j, k, l);
  status_x = _status[index];
  status_y = _status[index+_x*_y*_z*_t];
  status_z = _status[index+2*_x*_y*_z*_t];
}

irtkPoint irtkFreeFormTransformation4D::ControlPointLocation(int index) const
{
  int i, j, k, l;

  if (index >= _x*_y*_z*_t) {
    index -= _x*_y*_z*_t;
    if (index >= _x*_y*_z*_t) {
      index -= _x*_y*_z*_t;
    }
  }
  l = index/(_z*_y*_x);
  k = index%(_z*_y*_x)/(_y*_x);
  j = index%(_z*_y*_x)%(_y*_x)/_x;
  i = index%(_z*_y*_x)%(_y*_x)%_x;
  irtkPoint  p(i, j, k);
  this->LatticeToWorld(p);
  return p;
}

void irtkFreeFormTransformation4D::ControlPointLocation(int index, double &i, double &j, double &k) const
{
  int l;

  if (index >= _x*_y*_z*_t) {
    index -= _x*_y*_z*_t;
    if (index >= _x*_y*_z*_t) {
      index -= _x*_y*_z*_t;
    }
  }
  l = index/(_z*_y*_x);
  k = index%(_z*_y*_x)/(_y*_x);
  j = index%(_z*_y*_x)%(_y*_x)/_x;
  i = index%(_z*_y*_x)%(_y*_x)%_x;
  this->LatticeToWorld(i, j, k);
}

void irtkFreeFormTransformation4D::BoundingBox(irtkPoint &p1, irtkPoint &p2) const
{
  p1._x = 0;
  p1._y = 0;
  p1._z = 0;
  p2._x = _x - 1;
  p2._y = _y - 1;
  p2._z = _z - 1;
  this->LatticeToWorld(p1);
  this->LatticeToWorld(p2);
}

double irtkFreeFormTransformation4D::Inverse(double &x, double &y, double &z, double, double tolerance)
{
  int check;
  double error;

  // Initialize global variables
  irtkTransformationPointer  = this;
  x_invert = x;
  y_invert = y;
  z_invert = z;

  // Pointer to B-spline wrapper
  void (*Newton_function)(int, float [], float []) = irtkTransformationEvaluate;

  // Inverse
  float invert[4], f_invert[3];
  invert[0] = x;
  invert[1] = y;
  invert[2] = z;

  // Numerically approximate the inverse transformation
  newt2(invert-1, 3, &check, Newton_function);

  // Calculate error
  irtkTransformationEvaluate(3, invert-1, f_invert-1);
  error = sqrt(f_invert[0]*f_invert[0]+f_invert[1]*f_invert[1]+f_invert[2]*f_invert[2]);
  if (error > tolerance) {
    cout << "irtkFreeFormTransformation4D::Inverse: RMS error = " << error << "\n";
  }

  // Set output to solution
  x = invert[0];
  y = invert[1];
  z = invert[2];

  return error;
}
