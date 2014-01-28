/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#include <irtkTransformation.h>


irtkBSplineFreeFormTransformation4D::irtkBSplineFreeFormTransformation4D()
{
  int i;

  // Initialize control point domain
  _origin._x = 0.5;
  _origin._y = 0.5;
  _origin._z = 0.5;

  // Initialize control point domain (time)
  _tMin = 0;
  _tMax = 1;

  // Initialize x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Initialize y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Update z-axis
  _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
  _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
  _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];

  // Initialize control point dimensions
  _x = 2;
  _y = 2;
  _z = 2;
  _t = 2;

  // Initialize control point spacing
  _dx = 1;
  _dy = 1;
  _dz = 1;
  _dt = 1;

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Intialize memory for control point values
  _xdata = this->Allocate(_xdata, _x, _y, _z, _t);
  _ydata = this->Allocate(_ydata, _x, _y, _z, _t);
  _zdata = this->Allocate(_zdata, _x, _y, _z, _t);

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z*_t];
  for (i = 0; i < 3*_x*_y*_z*_t; i++) {
    _status[i] = _Active;
  }

  // Initialize lookup table
  _bspline.Initialize();
}

irtkBSplineFreeFormTransformation4D::irtkBSplineFreeFormTransformation4D(irtkBaseImage &image, double dx, double dy, double dz, double dt)
{
  int i;
  double x1, y1, z1, x2, y2, z2;

  // Figure out FOV
  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = image.GetX()-1;
  y2 = image.GetY()-1;
  z2 = image.GetZ()-1;
  image.ImageToWorld(x1, y1, z1);
  image.ImageToWorld(x2, y2, z2);

  // Initialize control point domain
  _origin._x = (x2 + x1) / 2.0;
  _origin._y = (y2 + y1) / 2.0;
  _origin._z = (z2 + z1) / 2.0;

  // Initialize control point domain
  _tMin = image.ImageToTime(0);
  _tMax = image.ImageToTime(image.GetT()-1);

  // Initialize x-axis and y-axis
  image.GetOrientation(_xaxis, _yaxis, _zaxis);

  double a = x1 * _xaxis[0] + y1 * _xaxis[1] + z1 * _xaxis[2];
  double b = x1 * _yaxis[0] + y1 * _yaxis[1] + z1 * _yaxis[2];
  double c = x1 * _zaxis[0] + y1 * _zaxis[1] + z1 * _zaxis[2];
  x1 = a;
  y1 = b;
  z1 = c;
  a = x2 * _xaxis[0] + y2 * _xaxis[1] + z2 * _xaxis[2];
  b = x2 * _yaxis[0] + y2 * _yaxis[1] + z2 * _yaxis[2];
  c = x2 * _zaxis[0] + y2 * _zaxis[1] + z2 * _zaxis[2];
  x2 = a;
  y2 = b;
  z2 = c;

  // Initialize control point dimensions
  _x = round((x2 - x1) / dx) + 1;
  _y = round((y2 - y1) / dy) + 1;
  _z = round((z2 - z1) / dz) + 1;
  _t = round((_tMax - _tMin) / dt) + 1;

  // Initialize control point spacing
  _dx = (_x > 1) ? ((x2 - x1)       / (_x - 1)) : dx;
  _dy = (_y > 1) ? ((y2 - y1)       / (_y - 1)) : dy;
  _dz = (_z > 1) ? ((z2 - z1)       / (_z - 1)) : dz;
  _dt = (_t > 1) ? ((_tMax - _tMin) / (_t - 1)) : dt;

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Intialize memory for control point values
  _xdata = this->Allocate(_xdata, _x, _y, _z, _t);
  _ydata = this->Allocate(_ydata, _x, _y, _z, _t);
  _zdata = this->Allocate(_zdata, _x, _y, _z, _t);

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z*_t];
  for (i = 0; i < 3*_x*_y*_z*_t; i++) {
    _status[i] = _Active;
  }

  // Initialize lookup table
  _bspline.Initialize();
}

irtkBSplineFreeFormTransformation4D::irtkBSplineFreeFormTransformation4D(irtkImageAttributes & attr, double dx, double dy, double dz, double dt)
{
  int i;
  double a, b, c, x1, y1, z1, x2, y2, z2;

  irtkMatrix matI2W = irtkBaseImage::GetImageToWorldMatrix(attr);
  irtkMatrix matW2I = irtkBaseImage::GetWorldToImageMatrix(attr);

  // Figure out FOV
  x1 = matI2W(0, 3);
  y1 = matI2W(1, 3);
  z1 = matI2W(2, 3);
  x2 = matI2W(0, 0)*(attr._x-1)+matI2W(0, 1)*(attr._y-1)+matI2W(0, 2)*(attr._z-1)+matI2W(0, 3);
  y2 = matI2W(1, 0)*(attr._x-1)+matI2W(1, 1)*(attr._y-1)+matI2W(1, 2)*(attr._z-1)+matI2W(1, 3);
  z2 = matI2W(2, 0)*(attr._x-1)+matI2W(2, 1)*(attr._y-1)+matI2W(2, 2)*(attr._z-1)+matI2W(2, 3);
  
  // Initialize control point domain
  _origin._x = (x2 + x1) / 2.0;
  _origin._y = (y2 + y1) / 2.0;
  _origin._z = (z2 + z1) / 2.0;
  
  // Initialize control point domain
  _tMin = attr._torigin;
  _tMax = attr._torigin + attr._t * attr._dt;
  
  // Initialize x-axis and y-axis
  _xaxis[0] = attr._xaxis[0];
  _xaxis[1] = attr._xaxis[1];
  _xaxis[2] = attr._xaxis[2];
  _yaxis[0] = attr._yaxis[0];
  _yaxis[1] = attr._yaxis[1];
  _yaxis[2] = attr._yaxis[2];
  _zaxis[0] = attr._zaxis[0];
  _zaxis[1] = attr._zaxis[1];
  _zaxis[2] = attr._zaxis[2];
  
  a = x1 * _xaxis[0] + y1 * _xaxis[1] + z1 * _xaxis[2];
  b = x1 * _yaxis[0] + y1 * _yaxis[1] + z1 * _yaxis[2];
  c = x1 * _zaxis[0] + y1 * _zaxis[1] + z1 * _zaxis[2];
  x1 = a;
  y1 = b;
  z1 = c;
  a = x2 * _xaxis[0] + y2 * _xaxis[1] + z2 * _xaxis[2];
  b = x2 * _yaxis[0] + y2 * _yaxis[1] + z2 * _yaxis[2];
  c = x2 * _zaxis[0] + y2 * _zaxis[1] + z2 * _zaxis[2];
  x2 = a;
  y2 = b;
  z2 = c;
  
  // Initialize control point dimensions
  _x = round((x2 - x1) / dx) + 1;
  _y = round((y2 - y1) / dy) + 1;
  _z = round((z2 - z1) / dz) + 1;
  _t = round((_tMax - _tMin) / dt) + 1;
  
  // Initialize control point spacing
  _dx = (_x > 1) ? ((x2 - x1)       / (_x - 1)) : dx;
  _dy = (_y > 1) ? ((y2 - y1)       / (_y - 1)) : dy;
  _dz = (_z > 1) ? ((z2 - z1)       / (_z - 1)) : dz;
  _dt = (_t > 1) ? ((_tMax - _tMin) / (_t - 1)) : dt;
  
  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);
  
  // Update transformation matrix
  this->UpdateMatrix();
  
  // Initialize memory for control point values
  _xdata = this->Allocate(_xdata, _x, _y, _z, _t);
  _ydata = this->Allocate(_ydata, _x, _y, _z, _t);
  _zdata = this->Allocate(_zdata, _x, _y, _z, _t);
  
  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z*_t];
  for (i = 0; i < 3*_x*_y*_z*_t; i++) {
    _status[i] = _Active;
  }
  
  // Initialize lookup table
  _bspline.Initialize();
}

irtkBSplineFreeFormTransformation4D::irtkBSplineFreeFormTransformation4D(double x1, double y1, double z1, double t1,
    double x2, double y2, double z2, double t2,
    double dx, double dy, double dz, double dt,
    double* xaxis, double* yaxis, double* zaxis)
{
  int i;

  // Initialize control point domain
  _origin._x = (x2 + x1) / 2.0;
  _origin._y = (y2 + y1) / 2.0;
  _origin._z = (z2 + z1) / 2.0;

  _tMin = t1;
  _tMax = t2;

  _xaxis[0] = xaxis[0];
  _xaxis[1] = xaxis[1];
  _xaxis[2] = xaxis[2];
  _yaxis[0] = yaxis[0];
  _yaxis[1] = yaxis[1];
  _yaxis[2] = yaxis[2];
  _zaxis[0] = zaxis[0];
  _zaxis[1] = zaxis[1];
  _zaxis[2] = zaxis[2];

  double a = x1 * _xaxis[0] + y1 * _xaxis[1] + z1 * _xaxis[2];
  double b = x1 * _yaxis[0] + y1 * _yaxis[1] + z1 * _yaxis[2];
  double c = x1 * _zaxis[0] + y1 * _zaxis[1] + z1 * _zaxis[2];
  x1 = a;
  y1 = b;
  z1 = c;
  a = x2 * _xaxis[0] + y2 * _xaxis[1] + z2 * _xaxis[2];
  b = x2 * _yaxis[0] + y2 * _yaxis[1] + z2 * _yaxis[2];
  c = x2 * _zaxis[0] + y2 * _zaxis[1] + z2 * _zaxis[2];
  x2 = a;
  y2 = b;
  z2 = c;

  // Initialize control point dimensions
  _x = round((x2 - x1) / dx) + 1;
  _y = round((y2 - y1) / dy) + 1;
  _z = round((z2 - z1) / dz) + 1;
  _t = round((_tMax - _tMin) / dt) + 1;

  // Initialize control point spacing
  _dx = (_x > 1) ? ((x2 - x1)       / (_x - 1)) : dx;
  _dy = (_y > 1) ? ((y2 - y1)       / (_y - 1)) : dy;
  _dz = (_z > 1) ? ((z2 - z1)       / (_z - 1)) : dz;
  _dt = (_t > 1) ? ((_tMax - _tMin) / (_t - 1)) : dt;

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Initialize memory for control point values
  _xdata = this->Allocate(_xdata, _x, _y, _z, _t);
  _ydata = this->Allocate(_ydata, _x, _y, _z, _t);
  _zdata = this->Allocate(_zdata, _x, _y, _z, _t);

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z*_t];
  for (i = 0; i < 3*_x*_y*_z*_t; i++) {
    _status[i] = _Active;
  }

  // Initialize lookup table
  _bspline.Initialize();
}

irtkBSplineFreeFormTransformation4D::irtkBSplineFreeFormTransformation4D(const irtkBSplineFreeFormTransformation4D &ffd) : irtkFreeFormTransformation4D(ffd)
{
  int i, j, k, l;

  // Initialize origin
  _origin = ffd._origin;

  // Initialize control point dimensions
  _x = ffd._x;
  _y = ffd._y;
  _z = ffd._z;
  _t = ffd._t;

  // Initialize control point spacing
  _dx = ffd._dx;
  _dy = ffd._dy;
  _dz = ffd._dz;
  _dt = ffd._dt;

  // Initialize time domain
  _tMin = ffd._tMin;
  _tMax = ffd._tMax;

  // Initialize x-axis
  _xaxis[0] = ffd._xaxis[0];
  _xaxis[1] = ffd._xaxis[1];
  _xaxis[2] = ffd._xaxis[2];

  // Initialize y-axis
  _yaxis[0] = ffd._yaxis[0];
  _yaxis[1] = ffd._yaxis[1];
  _yaxis[2] = ffd._yaxis[2];

  // Initialize z-axis
  _zaxis[0] = ffd._zaxis[0];
  _zaxis[1] = ffd._zaxis[1];
  _zaxis[2] = ffd._zaxis[2];

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Initialize memory for control point values
  _xdata = this->Allocate(_xdata, _x, _y, _z, _t);
  _ydata = this->Allocate(_ydata, _x, _y, _z, _t);
  _zdata = this->Allocate(_zdata, _x, _y, _z, _t);
  for (i = -4; i < _x+4; i++) {
    for (j = -4; j < _y+4; j++) {
      for (k = -4; k < _z+4; k++) {
        for (l = -4; l < _t+4; l++) {
          _xdata[l][k][j][i] = ffd._xdata[l][k][j][i];
          _ydata[l][k][j][i] = ffd._ydata[l][k][j][i];
          _zdata[l][k][j][i] = ffd._zdata[l][k][j][i];
        }
      }
    }
  }

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z*_t];
  for (i = 0; i < 3*_x*_y*_z*_t; i++) {
    _status[i] = ffd._status[i];
  }
}

irtkBSplineFreeFormTransformation4D::~irtkBSplineFreeFormTransformation4D()
{
  if (_xdata != NULL) _xdata = this->Deallocate(_xdata, _x, _y, _z, _t);
  if (_ydata != NULL) _ydata = this->Deallocate(_ydata, _x, _y, _z, _t);
  if (_zdata != NULL) _zdata = this->Deallocate(_zdata, _x, _y, _z, _t);
}

void irtkBSplineFreeFormTransformation4D::FFD2(double &x, double &y, double &z, double time) const
{
  double s, t, u, v, B_I, B_J, B_K, B_L;
  int a, b, c, d, i, j, k, l, S, T, U, V;

  // Now calculate the real stuff
  a = (int)floor(x)-1;
  b = (int)floor(y)-1;
  c = (int)floor(z)-1;
  d = (int)floor(time)-1;

  s = x-(a+1);
  t = y-(b+1);
  u = z-(c+1);
  v = time-(d+1);

  S = _bspline.VariableToIndex(s);
  T = _bspline.VariableToIndex(t);
  U = _bspline.VariableToIndex(u);
  V = _bspline.VariableToIndex(v);

  // Initialize displacement
  x = 0;
  y = 0;
  z = 0;

  for (l = 0; l < 4; l++) {
    B_L = _bspline.LookupTable[V][l];
    for (k = 0; k < 4; k++) {
      B_K = _bspline.LookupTable[U][k] * B_L;
      for (j = 0; j < 4; j++) {
        B_J = _bspline.LookupTable[T][j] * B_K;
        for (i = 0; i < 4; i++) {
          B_I = _bspline.LookupTable[S][i] * B_J;
          x += B_I * _xdata[d+l][c+k][b+j][a+i];
          y += B_I * _ydata[d+l][c+k][b+j][a+i];
          z += B_I * _zdata[d+l][c+k][b+j][a+i];
        }
      }
    }
  }
}

double irtkBSplineFreeFormTransformation4D::Approximate(double *x1, double *y1, double *z1, double *t1, double *x2, double *y2, double *z2, int no)
{
  int a, b, c, d, i, j, k, l, I, J, K, L, S, T, U, V, index;
  double s, t, u, v, x, y, z, B_I, B_J, B_K, B_L, basis, basis2, error, phi, norm, time;

  // Allocate memory
  double ****dx = NULL;
  double ****dy = NULL;
  double ****dz = NULL;
  double ****ds = NULL;
  dx = Allocate(dx, _x, _y, _z, _t);
  dy = Allocate(dy, _x, _y, _z, _t);
  dz = Allocate(dz, _x, _y, _z, _t);
  ds = Allocate(ds, _x, _y, _z, _t);

  // Subtract displacements which are approximated by current control points
  for (index = 0; index < no; index++) {
    x = x1[index];
    y = y1[index];
    z = z1[index];
    this->LocalDisplacement(x, y, z, t1[index]);
    x2[index] -= x;
    y2[index] -= y;
    z2[index] -= z;
  }

  // Initialize data structures
  for (l = -2; l < _t+2; l++) {
    for (k = -2; k < _z+2; k++) {
      for (j = -2; j < _y+2; j++) {
        for (i = -2; i < _x+2; i++) {
          dx[l][k][j][i] = 0;
          dy[l][k][j][i] = 0;
          dz[l][k][j][i] = 0;
          ds[l][k][j][i] = 0;
        }
      }
    }
  }

  // Initial loop: Calculate change of control points
  for (index = 0; index < no; index++) {
    x = x1[index];
    y = y1[index];
    z = z1[index];
    time = t1[index];
    this->WorldToLattice(x, y, z);
    time = this->TimeToLattice(time);
    a = (int)floor(x);
    b = (int)floor(y);
    c = (int)floor(z);
    d = (int)floor(time);
    s = x-a;
    t = y-b;
    u = z-c;
    v = time-d;
    S = _bspline.VariableToIndex(s);
    T = _bspline.VariableToIndex(t);
    U = _bspline.VariableToIndex(u);
    V = _bspline.VariableToIndex(v);
    norm = 0;
    for (l = 0; l < 4; l++) {
      B_L = _bspline.LookupTable[V][l];
      for (k = 0; k < 4; k++) {
        B_K = B_L * _bspline.LookupTable[U][k];
        for (j = 0; j < 4; j++) {
          B_J = B_K * _bspline.LookupTable[T][j];
          for (i = 0; i < 4; i++) {
            B_I = B_J * _bspline.LookupTable[S][i];
            norm += B_I * B_I;
          }
        }
      }
    }
    for (l = 0; l < 4; l++) {
      B_L = _bspline.LookupTable[V][l];
      L = l + d - 1;
      if ((L >= -2) && (L < _t+2)) {
        for (k = 0; k < 4; k++) {
          B_K = B_L * _bspline.LookupTable[U][k];
          K = k + c - 1;
          if ((K >= 0) && (K < _z+2)) {
            for (j = 0; j < 4; j++) {
              B_J = B_K * _bspline.LookupTable[T][j];
              J = j + b - 1;
              if ((J >= -2) && (J < _y+2)) {
                for (i = 0; i < 4; i++) {
                  B_I = B_J * _bspline.LookupTable[S][i];
                  I = i + a - 1;
                  if ((I >= -2) && (I < _x+2)) {
                    basis = B_I / norm;
                    basis2 = B_I * B_I;
                    phi = x2[index] * basis;
                    dx[L][K][J][I] += basis2 * phi;
                    phi = y2[index] * basis;
                    dy[L][K][J][I] += basis2 * phi;
                    phi = z2[index] * basis;
                    dz[L][K][J][I] += basis2 * phi;
                    ds[L][K][J][I] += basis2;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Add displacements which are approximated by current control points
  for (index = 0; index < no; index++) {
    x = x1[index];
    y = y1[index];
    z = z1[index];
    this->LocalDisplacement(x, y, z, t1[index]);
    x2[index] += x;
    y2[index] += y;
    z2[index] += z;
  }

  // Final loop: Calculate new control points
  for (l = -2; l < _t+2; l++) {
    for (k = -2; k < _z+2; k++) {
      for (j = -2; j < _y+2; j++) {
        for (i = -2; i < _x+2; i++) {
          if (ds[l][k][j][i] > 0) {
            _xdata[l][k][j][i] += dx[l][k][j][i] / ds[l][k][j][i];
            _ydata[l][k][j][i] += dy[l][k][j][i] / ds[l][k][j][i];
            _zdata[l][k][j][i] += dz[l][k][j][i] / ds[l][k][j][i];
          }
        }
      }
    }
  }

  // Calculate residual error
  error = 0;
  double max = 0;
  for (index = 0; index < no; index++) {
    x = x1[index];
    y = y1[index];
    z = z1[index];
    this->LocalDisplacement(x, y, z, t1[index]);
    x2[index] -= x;
    y2[index] -= y;
    z2[index] -= z;
    // Calculate error
    error += sqrt(x2[index]*x2[index]+y2[index]*y2[index]+z2[index]*z2[index]);
    if (sqrt(x2[index]*x2[index]+y2[index]*y2[index]+z2[index]*z2[index]) > max) max = sqrt(x2[index]*x2[index]+y2[index]*y2[index]+z2[index]*z2[index]);
  }
  error = error / (double)no;
  cout << max << endl;

  // Deallocate memory
  Deallocate(dx, _x, _y, _z, _t);
  Deallocate(dy, _x, _y, _z, _t);
  Deallocate(dz, _x, _y, _z, _t);
  Deallocate(ds, _x, _y, _z, _t);

  // Return error
  return error;
}

void irtkBSplineFreeFormTransformation4D::Jacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
  this->LocalJacobian(jac, x, y, z, t);
}

void irtkBSplineFreeFormTransformation4D::GlobalJacobian(irtkMatrix &jac, double, double, double, double)
{
  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Set matrix to identity
  jac(0, 0) = 1;
  jac(1, 1) = 1;
  jac(2, 2) = 1;
}

void irtkBSplineFreeFormTransformation4D::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
  int i, j, k, l;
  int floor_x, floor_y, floor_z, floor_t;
  int I, J, K, L;
  int IND_X, IND_Y, IND_Z, IND_T;

  double frac_x, frac_y, frac_z, frac_t;
  double coeff;
  double B_L, B_K, B_J, B_I, B_K_I, B_J_I, B_I_I;

  // The transformation maps (x, y, z, t) to (Tx, Ty, Tz)
  // Find the partial derivatives of the transformation Tx, Ty and Tz w.r.t x, y and z
  //     dTz/dz dTy/dz dTx/dz dTz/dy dTy/dy dTx/dy dTz/dx dTy/dx dTx/dx
  double z_k=0, y_k=0, x_k=0, z_j=0, y_j=0, x_j=0, z_i=0, y_i=0, x_i=0;

  // Partial derivatives with respect to time are not provided by this function.
  // I.e. none of dTx/dt , dTy/dt or dTz/dt is provided (they should be zero, right?).

  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  t = this->TimeToLattice(t);

  // Compute derivatives
  floor_x = (int)floor(x);
  floor_y = (int)floor(y);
  floor_z = (int)floor(z);
  floor_t = (int)floor(t);

  frac_x = x-floor_x;
  frac_y = y-floor_y;
  frac_z = z-floor_z;
  frac_t = t-floor_t;

  IND_X = _bspline.VariableToIndex(frac_x);
  IND_Y = _bspline.VariableToIndex(frac_y);
  IND_Z = _bspline.VariableToIndex(frac_z);
  IND_T = _bspline.VariableToIndex(frac_t);

  for (l = 0; l < 4; l++){
  	L = l + floor_t - 1;

  	if ((L >= 0) && (L < _t)){
			B_L   = _bspline.LookupTable[IND_T][l];
			// We are only returning the first three columns of the
			// Jacobian so do not need the following commented out bit
			// (the full Jacobian is a 4x3 matrix and we return a 3x3 one)

			// B_L_I = _bspline.LookupTable_I[IND_T][l];

  		for (k = 0; k < 4; k++) {
  			K = k + floor_z - 1;
  			if ((K >= 0) && (K < _z)) {
  				B_K   = _bspline.LookupTable[IND_Z][k];
  				B_K_I = _bspline.LookupTable_I[IND_Z][k];
  				for (j = 0; j < 4; j++) {
  					J = j + floor_y - 1;
  					if ((J >= 0) && (J < _y)) {
  						B_J   = _bspline.LookupTable[IND_Y][j];
  						B_J_I = _bspline.LookupTable_I[IND_Y][j];
  						for (i = 0; i < 4; i++) {
  							I = i + floor_x - 1;
  							if ((I >= 0) && (I < _x)) {
  								B_I   = _bspline.LookupTable[IND_X][i];
  								B_I_I = _bspline.LookupTable_I[IND_X][i];
  								coeff = B_I_I * B_J * B_K * B_L;
  								x_i += _xdata[L][K][J][I] * coeff;
  								y_i += _ydata[L][K][J][I] * coeff;
  								z_i += _zdata[L][K][J][I] * coeff;
  								coeff = B_I * B_J_I * B_K * B_L;
  								x_j += _xdata[L][K][J][I] * coeff;
  								y_j += _ydata[L][K][J][I] * coeff;
  								z_j += _zdata[L][K][J][I] * coeff;
  								coeff = B_I * B_J * B_K_I * B_L;
  								x_k += _xdata[L][K][J][I] * coeff;
  								y_k += _ydata[L][K][J][I] * coeff;
  								z_k += _zdata[L][K][J][I] * coeff;
  								// We are only returning the first three columns of the
  								// Jacobian so do not need the following
  								//  coeff = B_I * B_J * B_K * B_L_I;
  								// 	x_l += _xdata[L][K][J][I] * coeff;
  								//  y_l += _ydata[L][K][J][I] * coeff;
  								//  z_l += _zdata[L][K][J][I] * coeff;
  							}
  						}
  					}
  				}
  			}
  		}

  	}
  }

  // First 3 columns of the Jacobian matrix form a 3 x 3 sub-matrix
  jac.Initialize(3, 3);

  // Get deformation derivatives
  jac(0, 0) = x_i;
  jac(1, 0) = y_i;
  jac(2, 0) = z_i;
  jac(0, 1) = x_j;
  jac(1, 1) = y_j;
  jac(2, 1) = z_j;
  jac(0, 2) = x_k;
  jac(1, 2) = y_k;
  jac(2, 2) = z_k;

  // Convert derivatives to world coordinates
  jac = jac * _matW2L(0, 0, 3, 3);
  jac(0, 0) += 1;
  jac(1, 1) += 1;
  jac(2, 2) += 1;
}

void irtkBSplineFreeFormTransformation4D::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double t)
{
  // Transform index to control point location
  int i, j, k, l;
  this->IndexToLattice(dof, i, j, k, l);
  
  // Compute lattice coordinates
  this->WorldToLattice(x, y, z);
  t = this->TimeToLattice(t);
  
  // Calculate derivatives w.r.t. control point coordinates
  jac[0] = jac[1] = jac[2] = _bspline.B(x - i) * _bspline.B(y - j) * _bspline.B(z - k) * _bspline.B(t - l);
}

double irtkBSplineFreeFormTransformation4D::Bending(double, double, double, double)
{
  cerr << "irtkBSplineFreeFormTransformation4D::Bending: Not implemented yet" << endl;
  exit(1);
}

int irtkBSplineFreeFormTransformation4D::CheckHeader(char *)
{
  return false;
}

void irtkBSplineFreeFormTransformation4D::Interpolate(double *, double *, double *)
{
  cerr << "irtkBSplineFreeFormTransformation4D::Interpolate: Not implemented yet" << endl;
  exit(1);
}

void irtkBSplineFreeFormTransformation4D::Subdivide2D()
{
  int i, j, k, l, i1, j1, i2, j2;

  // Weights for subdivision
  double w[2][3];
  w[1][0] = 0;
  w[1][1] = 1.0/2.0;
  w[1][2] = 1.0/2.0;
  w[0][0] = 1.0/8.0;
  w[0][1] = 6.0/8.0;
  w[0][2] = 1.0/8.0;

  // Allocate memory for new control points
  double ****x = NULL;
  double ****y = NULL;
  double ****z = NULL;
  x = this->Allocate(x, 2*_x-1, 2*_y-1, _z, _t);
  y = this->Allocate(y, 2*_x-1, 2*_y-1, _z, _t);
  z = this->Allocate(z, 2*_x-1, 2*_y-1, _z, _t);

  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        for (l = 0; l < _t; l++) {
          for (i1 = 0; i1 < 2; i1++) {
            for (j1 = 0; j1 < 2; j1++) {
              x[l][k][2*j+j1][2*i+i1] = 0;
              y[l][k][2*j+j1][2*i+i1] = 0;
              z[l][k][2*j+j1][2*i+i1] = 0;
              for (i2 = 0; i2 < 3; i2++) {
                for (j2 = 0; j2 < 3; j2++) {
                  x[l][k][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * _xdata[l][k][j+j2-1][i+i2-1];
                  y[l][k][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * _ydata[l][k][j+j2-1][i+i2-1];
                  z[l][k][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * _zdata[l][k][j+j2-1][i+i2-1];
                }
              }
            }
          }
        }
      }
    }
  }

  // Deallocate points
  this->Deallocate(_xdata, _x, _y, _z, _t);
  this->Deallocate(_ydata, _x, _y, _z, _t);
  this->Deallocate(_zdata, _x, _y, _z, _t);
  delete []_status;

  // Update pointers to control points
  _xdata = x;
  _ydata = y;
  _zdata = z;

  // Increase number of control points
  _x = 2*_x - 1;
  _y = 2*_y - 1;

  // Recalculate control point spacing
  _dx /= 2.0;
  _dy /= 2.0;

  // Update transformation matrix
  this->UpdateMatrix();

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z*_t];
  for (i = 0; i < 3*_x*_y*_z*_t; i++) {
    _status[i] = _Active;
  }
}

void irtkBSplineFreeFormTransformation4D::Subdivide3D()
{
  int i, j, k, l, i1, j1, k1, i2, j2, k2;

  // Weights for subdivision
  double w[2][3];
  w[1][0] = 0;
  w[1][1] = 1.0/2.0;
  w[1][2] = 1.0/2.0;
  w[0][0] = 1.0/8.0;
  w[0][1] = 6.0/8.0;
  w[0][2] = 1.0/8.0;

  // Allocate memory for new control points
  double ****x = NULL;
  double ****y = NULL;
  double ****z = NULL;
  x = this->Allocate(x, 2*_x-1, 2*_y-1, 2*_z-1, _t);
  y = this->Allocate(y, 2*_x-1, 2*_y-1, 2*_z-1, _t);
  z = this->Allocate(z, 2*_x-1, 2*_y-1, 2*_z-1, _t);

  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        for (l = 0; l < _t; l++) {
          for (i1 = 0; i1 < 2; i1++) {
            for (j1 = 0; j1 < 2; j1++) {
              for (k1 = 0; k1 < 2; k1++) {
                x[l][2*k+k1][2*j+j1][2*i+i1] = 0;
                y[l][2*k+k1][2*j+j1][2*i+i1] = 0;
                z[l][2*k+k1][2*j+j1][2*i+i1] = 0;
                for (i2 = 0; i2 < 3; i2++) {
                  for (j2 = 0; j2 < 3; j2++) {
                    for (k2 = 0; k2 < 3; k2++) {
                      x[l][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * _xdata[l][k+k2-1][j+j2-1][i+i2-1];
                      y[l][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * _ydata[l][k+k2-1][j+j2-1][i+i2-1];
                      z[l][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * _zdata[l][k+k2-1][j+j2-1][i+i2-1];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Deallocate points
  this->Deallocate(_xdata, _x, _y, _z, _t);
  this->Deallocate(_ydata, _x, _y, _z, _t);
  this->Deallocate(_zdata, _x, _y, _z, _t);
  delete []_status;

  // Update pointers to control points
  _xdata = x;
  _ydata = y;
  _zdata = z;

  // Increase number of control points
  _x = 2*_x - 1;
  _y = 2*_y - 1;
  _z = 2*_z - 1;

  // Recalculate control point spacing
  _dx /= 2.0;
  _dy /= 2.0;
  _dz /= 2.0;

  // Update transformation matrix
  this->UpdateMatrix();

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z*_t];
  for (i = 0; i < 3*_x*_y*_z*_t; i++) {
    _status[i] = _Active;
  }
}

void irtkBSplineFreeFormTransformation4D::Subdivide4D()
{
  int i, j, k, l, i1, j1, k1, l1, i2, j2, k2, l2;

  // Weights for subdivision
  double w[2][3];
  w[1][0] = 0;
  w[1][1] = 1.0/2.0;
  w[1][2] = 1.0/2.0;
  w[0][0] = 1.0/8.0;
  w[0][1] = 6.0/8.0;
  w[0][2] = 1.0/8.0;

  // Allocate memory for new control points
  double ****x = NULL;
  double ****y = NULL;
  double ****z = NULL;
  x = this->Allocate(x, 2*_x-1, 2*_y-1, 2*_z-1, 2*_t-1);
  y = this->Allocate(y, 2*_x-1, 2*_y-1, 2*_z-1, 2*_t-1);
  z = this->Allocate(z, 2*_x-1, 2*_y-1, 2*_z-1, 2*_t-1);

  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        for (l = 0; l < _t; l++) {
          for (i1 = 0; i1 < 2; i1++) {
            for (j1 = 0; j1 < 2; j1++) {
              for (k1 = 0; k1 < 2; k1++) {
                for (l1 = 0; l1 < 2; l1++) {
                  x[2*l+l1][2*k+k1][2*j+j1][2*i+i1] = 0;
                  y[2*l+l1][2*k+k1][2*j+j1][2*i+i1] = 0;
                  z[2*l+l1][2*k+k1][2*j+j1][2*i+i1] = 0;
                  for (i2 = 0; i2 < 3; i2++) {
                    for (j2 = 0; j2 < 3; j2++) {
                      for (k2 = 0; k2 < 3; k2++) {
                        for (l2 = 0; l2 < 3; l2++) {
                          x[2*l+l1][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * w[l1][l2] * _xdata[l+l2-1][k+k2-1][j+j2-1][i+i2-1];
                          y[2*l+l1][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * w[l1][l2] * _ydata[l+l2-1][k+k2-1][j+j2-1][i+i2-1];
                          z[2*l+l1][2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] * w[j1][j2] * w[k1][k2] * w[l1][l2] * _zdata[l+l2-1][k+k2-1][j+j2-1][i+i2-1];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Deallocate points
  this->Deallocate(_xdata, _x, _y, _z, _t);
  this->Deallocate(_ydata, _x, _y, _z, _t);
  this->Deallocate(_zdata, _x, _y, _z, _t);
  delete []_status;

  // Update pointers to control points
  _xdata = x;
  _ydata = y;
  _zdata = z;

  // Increase number of control points
  _x = 2*_x - 1;
  _y = 2*_y - 1;
  _z = 2*_z - 1;
  _t = 2*_t - 1;

  // Recalculate control point spacing
  _dx /= 2.0;
  _dy /= 2.0;
  _dz /= 2.0;
  _dt /= 2.0;

  // Update transformation matrix
  this->UpdateMatrix();

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z*_t];
  for (i = 0; i < 3*_x*_y*_z*_t; i++) {
    _status[i] = _Active;
  }
}

void irtkBSplineFreeFormTransformation4D::BoundingBoxCP(int index, irtkPoint &p1, irtkPoint &p2, double fraction) const
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
  p1 = irtkPoint(i-2*fraction, j-2*fraction, k-2*fraction);
  p2 = irtkPoint(i+2*fraction, j+2*fraction, k+2*fraction);
  this->LatticeToWorld(p1);
  this->LatticeToWorld(p2);
}

void irtkBSplineFreeFormTransformation4D::BoundingBoxCP(int index, double &x1, double &y1, double &z1, double &t1, double &x2, double &y2, double &z2, double &t2, double fraction) const
{
  if (index >= _x*_y*_z*_t) {
    index -= _x*_y*_z*_t;
    if (index >= _x*_y*_z*_t) {
      index -= _x*_y*_z*_t;
    }
  }
  x1 = index%(_z*_y*_x)%(_y*_x)%_x-2*fraction;
  y1 = index%(_z*_y*_x)%(_y*_x)/_x-2*fraction;
  z1 = index%(_z*_y*_x)/(_y*_x)-2*fraction;
  t1 = index/(_z*_y*_x)-2*fraction;
  x2 = index%(_z*_y*_x)%(_y*_x)%_x+2*fraction;
  y2 = index%(_z*_y*_x)%(_y*_x)/_x+2*fraction;
  z2 = index%(_z*_y*_x)/(_y*_x)+2*fraction;
  t2 = index/(_z*_y*_x)+2*fraction;
  this->LatticeToWorld(x1, y1, z1);
  this->LatticeToWorld(x2, y2, z2);
  t1 = this->LatticeToTime(t1);
  t2 = this->LatticeToTime(t2);
}

void irtkBSplineFreeFormTransformation4D::BoundingBoxImage(irtkGreyImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
{
	irtkPoint p1, p2;
  double x1, y1, z1, x2, y2, z2;

  // Calculate bounding box in world coordinates
  this->BoundingBoxCP(index, p1, p2, fraction);
  x1 = p1._x;
  y1 = p1._y;
  z1 = p1._z;
  x2 = p2._x;
  y2 = p2._y;
  z2 = p2._z;

  // Transform world coordinates to image coordinates
  image->WorldToImage(x1, y1, z1);
  image->WorldToImage(x2, y2, z2);

  // Calculate bounding box in image coordinates
  i1 = (x1 < 0) ? 0 : int(x1)+1;
  j1 = (y1 < 0) ? 0 : int(y1)+1;
  k1 = (z1 < 0) ? 0 : int(z1)+1;
  i2 = (int(x2) >= image->GetX()) ? image->GetX()-1 : int(x2);
  j2 = (int(y2) >= image->GetY()) ? image->GetY()-1 : int(y2);
  k2 = (int(z2) >= image->GetZ()) ? image->GetZ()-1 : int(z2);
}

void irtkBSplineFreeFormTransformation4D::BoundingBoxImage(irtkRealImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
{
	irtkPoint p1, p2;
  double x1, y1, z1, x2, y2, z2;

  // Calculate bounding box in world coordinates
  this->BoundingBoxCP(index, p1, p2, fraction);
  x1 = p1._x;
  y1 = p1._y;
  z1 = p1._z;
  x2 = p2._x;
  y2 = p2._y;
  z2 = p2._z;

  // Transform world coordinates to image coordinates
  image->WorldToImage(x1, y1, z1);
  image->WorldToImage(x2, y2, z2);

  // Calculate bounding box in image coordinates
  i1 = (x1 < 0) ? 0 : int(x1)+1;
  j1 = (y1 < 0) ? 0 : int(y1)+1;
  k1 = (z1 < 0) ? 0 : int(z1)+1;
  i2 = (int(x2) >= image->GetX()) ? image->GetX()-1 : int(x2);
  j2 = (int(y2) >= image->GetY()) ? image->GetY()-1 : int(y2);
  k2 = (int(z2) >= image->GetZ()) ? image->GetZ()-1 : int(z2);
}

void irtkBSplineFreeFormTransformation4D::Print()
{
  // Print keyword and no. of DOFs
  cout << "BSplineFFD4D: " << this->NumberOfDOFs() << endl;

  // Print no. of control points
  cout << "Control points: " << _x << " x " << _y << " x " << _z << " x " << _t << endl;
  cout << "Spacing: " << _dx << " x " << _dy << " x " << _dz << " x " << _dt << endl;
  cout << "Time: " << _tMin << " to " << _tMax << endl;
  cout << "Origin: " << _origin._x << " " << _origin._y << " " << _origin._z << endl;
  // Print x-axis
  cout << "X-axis is " << _xaxis[0] << " " << _xaxis[1] << " " << _xaxis[2] << endl;
  // Print x-axis
  cout << "Y-axis is " << _yaxis[0] << " " << _yaxis[1] << " " << _yaxis[2] << endl;
  // Print x-axis
  cout << "Z-axis is " << _zaxis[0] << " " << _zaxis[1] << " " << _zaxis[2] << endl;
}

irtkCifstream& irtkBSplineFreeFormTransformation4D::Read(irtkCifstream& from)
{
  double *data;
  int i, j, k, l, index;
  unsigned int magic_no, trans_type;

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkBSplineFreeFormTransformation4D::Read: Not a vaild transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if (trans_type != IRTKTRANSFORMATION_BSPLINE_FFD_4D) {
    cerr << "irtkBSplineFreeFormTransformation4D::Read: Not a vaild B-Spline FFD transformation" << endl;
    exit(1);
  }

  // Free memory if necessary
  _xdata = Deallocate(_xdata, _x, _y, _z, _t);
  _ydata = Deallocate(_ydata, _x, _y, _z, _t);
  _zdata = Deallocate(_zdata, _x, _y, _z, _t);
  delete []_status;

  // Read no of control points
  from.ReadAsInt(&_x, 1);
  from.ReadAsInt(&_y, 1);
  from.ReadAsInt(&_z, 1);
  from.ReadAsInt(&_t, 1);

  // Read orientation of bounding box
  from.ReadAsDouble(_xaxis, 3);
  from.ReadAsDouble(_yaxis, 3);
  from.ReadAsDouble(_zaxis, 3);

  // Read spacing of bounding box
  from.ReadAsDouble(&_dx, 1);
  from.ReadAsDouble(&_dy, 1);
  from.ReadAsDouble(&_dz, 1);
  from.ReadAsDouble(&_dt, 1);

  // Read spacing of bounding box
  from.ReadAsDouble(&_origin._x, 1);
  from.ReadAsDouble(&_origin._y, 1);
  from.ReadAsDouble(&_origin._z, 1);

  // Read time domain
  from.ReadAsDouble(&_tMin, 1);
  from.ReadAsDouble(&_tMax, 1);

  // Initialize control points
  _xdata = Allocate(_xdata, _x, _y, _z, _t);
  _ydata = Allocate(_ydata, _x, _y, _z, _t);
  _zdata = Allocate(_zdata, _x, _y, _z, _t);

  // Initialize control points
  for (i = -2; i < _x+2; i++) {
    for (j = -2; j < _y+2; j++) {
      for (k = -2; k < _z+2; k++) {
        for (l = -2; l < _t+2; l++) {
          _xdata[l][k][j][i] = 0;
          _ydata[l][k][j][i] = 0;
          _zdata[l][k][j][i] = 0;
        }
      }
    }
  }

  // Allocate temporary memory
  data = new double[3*_x*_y*_z*_t];

  // Read control point data
  from.ReadAsDouble(data, 3*_x*_y*_z*_t);

  // Convert data
  index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        for (l = 0; l < _t; l++) {
          _xdata[l][k][j][i] = data[index];
          _ydata[l][k][j][i] = data[index+1];
          _zdata[l][k][j][i] = data[index+2];
          index += 3;
        }
      }
    }
  }

  // Free temporary memory
  delete []data;

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z*_t];

  // Read control point status
  from.ReadAsInt((int *)_status, 3*_x*_y*_z*_t);

  // Update transformation matrix
  this->UpdateMatrix();

  return from;
}

irtkCofstream& irtkBSplineFreeFormTransformation4D::Write(irtkCofstream& to)
{
  double *data;
  int i, j, k, l, index;
  unsigned int magic_no, trans_type;

  // Write magic no. for transformations
  magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  trans_type = IRTKTRANSFORMATION_BSPLINE_FFD_4D;
  to.WriteAsUInt(&trans_type, 1);

  // Write no of control points
  to.WriteAsInt(&_x, 1);
  to.WriteAsInt(&_y, 1);
  to.WriteAsInt(&_z, 1);
  to.WriteAsInt(&_t, 1);

  // Write orientation of bounding box
  to.WriteAsDouble(_xaxis, 3);
  to.WriteAsDouble(_yaxis, 3);
  to.WriteAsDouble(_zaxis, 3);

  // Write spacing of bounding box
  to.WriteAsDouble(&_dx, 1);
  to.WriteAsDouble(&_dy, 1);
  to.WriteAsDouble(&_dz, 1);
  to.WriteAsDouble(&_dt, 1);

  // Write spacing of bounding box
  to.WriteAsDouble(&_origin._x, 1);
  to.WriteAsDouble(&_origin._y, 1);
  to.WriteAsDouble(&_origin._z, 1);

  // Write time domain
  to.WriteAsDouble(&_tMin, 1);
  to.WriteAsDouble(&_tMax, 1);

  // Allocate temporary memory
  data = new double[3*_x*_y*_z*_t];

  // Convert data
  index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        for (l = 0; l < _t; l++) {
          data[index]   = _xdata[l][k][j][i];
          data[index+1] = _ydata[l][k][j][i];
          data[index+2] = _zdata[l][k][j][i];
          index += 3;
        }
      }
    }
  }

  // Write control point data
  to.WriteAsDouble(data, 3*_x*_y*_z*_t);

  // Free temporary memory
  delete []data;

  // Write control point status
  to.WriteAsInt((int *)_status, 3*_x*_y*_z*_t);

  return to;
}
