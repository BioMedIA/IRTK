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

#define LUTSIZE (double)(FFDLOOKUPTABLESIZE-1)

double irtkBSplineFreeFormTransformation4D::LookupTable   [FFDLOOKUPTABLESIZE][4];

double irtkBSplineFreeFormTransformation4D::LookupTable_I [FFDLOOKUPTABLESIZE][4];

double irtkBSplineFreeFormTransformation4D::LookupTable_II[FFDLOOKUPTABLESIZE][4];

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
  for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
    this->LookupTable[i][0]   = this->B0(i/LUTSIZE);
    this->LookupTable[i][1]   = this->B1(i/LUTSIZE);
    this->LookupTable[i][2]   = this->B2(i/LUTSIZE);
    this->LookupTable[i][3]   = this->B3(i/LUTSIZE);
    this->LookupTable_I[i][0] = this->B0_I(i/LUTSIZE);
    this->LookupTable_I[i][1] = this->B1_I(i/LUTSIZE);
    this->LookupTable_I[i][2] = this->B2_I(i/LUTSIZE);
    this->LookupTable_I[i][3] = this->B3_I(i/LUTSIZE);
    this->LookupTable_II[i][0] = this->B0_II(i/LUTSIZE);
    this->LookupTable_II[i][1] = this->B1_II(i/LUTSIZE);
    this->LookupTable_II[i][2] = this->B2_II(i/LUTSIZE);
    this->LookupTable_II[i][3] = this->B3_II(i/LUTSIZE);
  }
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
  _dx = (x2 - x1) / (_x - 1);
  _dy = (y2 - y1) / (_y - 1);
  _dz = (z2 - z1) / (_z - 1);
  _dt = (_tMax - _tMin) / (_t - 1);

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
  for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
    this->LookupTable[i][0]   = this->B0(i/LUTSIZE);
    this->LookupTable[i][1]   = this->B1(i/LUTSIZE);
    this->LookupTable[i][2]   = this->B2(i/LUTSIZE);
    this->LookupTable[i][3]   = this->B3(i/LUTSIZE);
    this->LookupTable_I[i][0] = this->B0_I(i/LUTSIZE);
    this->LookupTable_I[i][1] = this->B1_I(i/LUTSIZE);
    this->LookupTable_I[i][2] = this->B2_I(i/LUTSIZE);
    this->LookupTable_I[i][3] = this->B3_I(i/LUTSIZE);
    this->LookupTable_II[i][0] = this->B0_II(i/LUTSIZE);
    this->LookupTable_II[i][1] = this->B1_II(i/LUTSIZE);
    this->LookupTable_II[i][2] = this->B2_II(i/LUTSIZE);
    this->LookupTable_II[i][3] = this->B3_II(i/LUTSIZE);
  }
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
  _dx = (x2 - x1) / (_x - 1);
  _dy = (y2 - y1) / (_y - 1);
  _dz = (z2 - z1) / (_z - 1);
  _dt = (_tMax - _tMin) / (_t - 1);

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
  for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
    this->LookupTable[i][0]   = this->B0(i/LUTSIZE);
    this->LookupTable[i][1]   = this->B1(i/LUTSIZE);
    this->LookupTable[i][2]   = this->B2(i/LUTSIZE);
    this->LookupTable[i][3]   = this->B3(i/LUTSIZE);
    this->LookupTable_I[i][0] = this->B0_I(i/LUTSIZE);
    this->LookupTable_I[i][1] = this->B1_I(i/LUTSIZE);
    this->LookupTable_I[i][2] = this->B2_I(i/LUTSIZE);
    this->LookupTable_I[i][3] = this->B3_I(i/LUTSIZE);
    this->LookupTable_II[i][0] = this->B0_II(i/LUTSIZE);
    this->LookupTable_II[i][1] = this->B1_II(i/LUTSIZE);
    this->LookupTable_II[i][2] = this->B2_II(i/LUTSIZE);
    this->LookupTable_II[i][3] = this->B3_II(i/LUTSIZE);
  }
}

irtkBSplineFreeFormTransformation4D::irtkBSplineFreeFormTransformation4D(const irtkBSplineFreeFormTransformation4D &ffd)
{
  int i, j, k, l;

  // Initialize origin
  _origin = ffd._origin;

  // Initialize control point dimensions
  _x = ffd._x;
  _y = ffd._y;
  _z = ffd._z;
  _t = ffd._t;

  // Intialize control point spacing
  _dx = ffd._dx;
  _dy = ffd._dy;
  _dz = ffd._dz;
  _dt = ffd._dt;

  // Intialize time domain
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

  // Initialize lookup table
  for (i = 0; i < FFDLOOKUPTABLESIZE; i++) {
    this->LookupTable[i][0]   = this->B0(i/LUTSIZE);
    this->LookupTable[i][1]   = this->B1(i/LUTSIZE);
    this->LookupTable[i][2]   = this->B2(i/LUTSIZE);
    this->LookupTable[i][3]   = this->B3(i/LUTSIZE);
    this->LookupTable_I[i][0] = this->B0_I(i/LUTSIZE);
    this->LookupTable_I[i][1] = this->B1_I(i/LUTSIZE);
    this->LookupTable_I[i][2] = this->B2_I(i/LUTSIZE);
    this->LookupTable_I[i][3] = this->B3_I(i/LUTSIZE);
    this->LookupTable_II[i][0] = this->B0_II(i/LUTSIZE);
    this->LookupTable_II[i][1] = this->B1_II(i/LUTSIZE);
    this->LookupTable_II[i][2] = this->B2_II(i/LUTSIZE);
    this->LookupTable_II[i][3] = this->B3_II(i/LUTSIZE);
  }
}

irtkBSplineFreeFormTransformation4D::~irtkBSplineFreeFormTransformation4D()
{
  // Free memory for control points if necessary
  if (_xdata != NULL) _xdata = this->Deallocate(_xdata, _x, _y, _z, _t);
  if (_ydata != NULL) _ydata = this->Deallocate(_ydata, _x, _y, _z, _t);
  if (_zdata != NULL) _zdata = this->Deallocate(_zdata, _x, _y, _z, _t);

  _x = 0;
  _y = 0;
  _z = 0;
}

void irtkBSplineFreeFormTransformation4D::FFD1(double &x, double &y, double &z, double time) const
{
  double s, t, u, v, B_I, B_J, B_K, B_L;
  int a, b, c, d, i, j, k, l, S, T, U, V;

  // Check if there is some work to do
  if ((x < -2) || (y < -2) || (z < -2) || (time < -2) || (x > _x+1) || (y > _y+1) || (z > _z+1) || (time > _t+1)) {
    x = 0;
    y = 0;
    z = 0;
    return;
  }

  // Now calculate the real stuff
  a = (int)floor(x)-1;
  b = (int)floor(y)-1;
  c = (int)floor(z)-1;
  d = (int)floor(time)-1;

  s = x-(a+1);
  t = y-(b+1);
  u = z-(c+1);
  v = time-(d+1);

  S = round(LUTSIZE*s);
  T = round(LUTSIZE*t);
  U = round(LUTSIZE*u);
  V = round(LUTSIZE*v);

  // Initialize displacement
  x = 0;
  y = 0;
  z = 0;

  for (l = 0; l < 4; l++) {
    B_L = this->LookupTable[V][l];
    for (k = 0; k < 4; k++) {
      B_K = this->LookupTable[U][k] * B_L;
      for (j = 0; j < 4; j++) {
        B_J = this->LookupTable[T][j] * B_K;
        for (i = 0; i < 4; i++) {
          B_I = this->LookupTable[S][i] * B_J;
          x += B_I * _xdata[d+l][c+k][b+j][a+i];
          y += B_I * _ydata[d+l][c+k][b+j][a+i];
          z += B_I * _zdata[d+l][c+k][b+j][a+i];
        }
      }
    }
  }
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

  S = round(LUTSIZE*s);
  T = round(LUTSIZE*t);
  U = round(LUTSIZE*u);
  V = round(LUTSIZE*v);

  // Initialize displacement
  x = 0;
  y = 0;
  z = 0;

  for (l = 0; l < 4; l++) {
    B_L = this->LookupTable[V][l];
    for (k = 0; k < 4; k++) {
      B_K = this->LookupTable[U][k] * B_L;
      for (j = 0; j < 4; j++) {
        B_J = this->LookupTable[T][j] * B_K;
        for (i = 0; i < 4; i++) {
          B_I = this->LookupTable[S][i] * B_J;
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
    S = round(LUTSIZE*s);
    T = round(LUTSIZE*t);
    U = round(LUTSIZE*u);
    V = round(LUTSIZE*v);
    norm = 0;
    for (l = 0; l < 4; l++) {
      B_L = this->LookupTable[V][l];
      for (k = 0; k < 4; k++) {
        B_K = B_L * this->LookupTable[U][k];
        for (j = 0; j < 4; j++) {
          B_J = B_K * this->LookupTable[T][j];
          for (i = 0; i < 4; i++) {
            B_I = B_J * this->LookupTable[S][i];
            norm += B_I * B_I;
          }
        }
      }
    }
    for (l = 0; l < 4; l++) {
      B_L = this->LookupTable[V][l];
      L = l + d - 1;
      if ((L >= -2) && (L < _t+2)) {
        for (k = 0; k < 4; k++) {
          B_K = B_L * this->LookupTable[U][k];
          K = k + c - 1;
          if ((K >= 0) && (K < _z+2)) {
            for (j = 0; j < 4; j++) {
              B_J = B_K * this->LookupTable[T][j];
              J = j + b - 1;
              if ((J >= -2) && (J < _y+2)) {
                for (i = 0; i < 4; i++) {
                  B_I = B_J * this->LookupTable[S][i];
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

void irtkBSplineFreeFormTransformation4D::GlobalJacobian(irtkMatrix &jac, double x, double y, double z, double t)
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
  cerr << "irtkBSplineFreeFormTransformation4D::LocalJacobian: Not implemented yet" << endl;
  exit(1);
}

double irtkBSplineFreeFormTransformation4D::Bending(double x, double y, double z, double t)
{
  cerr << "irtkBSplineFreeFormTransformation4D::Bending: Not implemented yet" << endl;
  exit(1);
}

int irtkBSplineFreeFormTransformation4D::CheckHeader(char *name)
{
  return False;
}

void irtkBSplineFreeFormTransformation4D::Interpolate(double *, double *, double *)
{
  cerr << "irtkBSplineFreeFormTransformation4D::Interpolate: Not implemented yet" << endl;
  exit(1);
}

void irtkBSplineFreeFormTransformation4D::Subdivide()
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

void irtkBSplineFreeFormTransformation4D::Subdivide2D()
{
  cerr << "irtkBSplineFreeFormTransformation4D::Subdivide2D: Not implemented yet" << endl;
  exit(1);
}

void irtkBSplineFreeFormTransformation4D::Subdivide3D()
{
  cerr << "irtkBSplineFreeFormTransformation4D::Subdivide3D: Not implemented yet" << endl;
  exit(1);
}

void irtkBSplineFreeFormTransformation4D::BoundingBox(int index, irtkPoint &p1, irtkPoint &p2, double fraction) const
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

void irtkBSplineFreeFormTransformation4D::BoundingBox(int index, double &x1, double &y1, double &z1, double &t1, double &x2, double &y2, double &z2, double &t2, double fraction) const
{
  if (index >= _x*_y*_z*_t) {
    index -= _x*_y*_z*_t;
    if (index >= _x*_y*_z*_t) {
      index -= _x*_y*_z*_t;
    }
  }
  x1 = index/(_z*_y*_x)-2*fraction;
  y1 = index%(_z*_y*_x)/(_y*_x)-2*fraction;
  z1 = index%(_z*_y*_x)%(_y*_x)/_x-2*fraction;
  t1 = index%(_z*_y*_x)%(_y*_x)%_x-2*fraction;
  x2 = index/(_z*_y*_x)+2*fraction;
  y2 = index%(_z*_y*_x)/(_y*_x)+2*fraction;
  z2 = index%(_z*_y*_x)%(_y*_x)/_x+2*fraction;
  t2 = index%(_z*_y*_x)%(_y*_x)%_x+2*fraction;
  this->LatticeToWorld(x1, y1, z1);
  this->LatticeToWorld(x2, y2, z2);
  t1 = this->LatticeToTime(t1);
  t2 = this->LatticeToTime(t2);
}

void irtkBSplineFreeFormTransformation4D::BoundingBox(irtkGreyImage *image, int index, int &i1, int &j1, int &k1, int &l1, int &i2, int &j2, int &k2, int &l2, double fraction) const
{
  double x1, y1, z1, t1, x2, y2, z2, t2;

  // Calculate bounding box in world coordinates
  this->BoundingBox(index, x1, y1, z1, t1, x2, y2, z2, t2, fraction);

  // Transform world coordinates to image coordinates
  image->WorldToImage(x1, y1, z1);
  image->WorldToImage(x2, y2, z2);
  t1 = image->TimeToImage(t1);
  t2 = image->TimeToImage(t1);

  // Calculate bounding box in image coordinates
  i1 = (x1 < 0) ? 0 : int(x1)+1;
  j1 = (y1 < 0) ? 0 : int(y1)+1;
  k1 = (z1 < 0) ? 0 : int(z1)+1;
  l1 = (t1 < 0) ? 0 : int(t1)+1;
  i2 = (int(x2) >= image->GetX()) ? image->GetX()-1 : int(x2);
  j2 = (int(y2) >= image->GetY()) ? image->GetY()-1 : int(y2);
  k2 = (int(z2) >= image->GetZ()) ? image->GetZ()-1 : int(z2);
  l2 = (int(t2) >= image->GetT()) ? image->GetT()-1 : int(t2);
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
