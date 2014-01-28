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

#include <newt2.h>

irtkLinearFreeFormTransformation::irtkLinearFreeFormTransformation()
{
  int i;

  // Initialize control point domain
  _origin._x = 0.5;
  _origin._y = 0.5;
  _origin._z = 0.5;

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

  // Initialize control point spacing
  _dx = 1;
  _dy = 1;
  _dz = 1;

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Intialize memory for control point values
  _data = this->Allocate(_data, _x, _y, _z);

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z];
  if (_z == 1) {
    for (i = 0; i < 2*_x*_y*_z; i++) {
      _status[i] = _Active;
    }
    for (i = 2*_x*_y*_z; i < 3*_x*_y*_z; i++) {
      _status[i] = _Passive;
    }
  } else {
    for (i = 0; i < 3*_x*_y*_z; i++) {
      _status[i] = _Active;
    }
  }
}

irtkLinearFreeFormTransformation::irtkLinearFreeFormTransformation(irtkBaseImage &image, double dx, double dy, double dz)
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
  if (z2 > z1) {
    _z = round((z2 - z1) / dz) + 1;
  } else {
    _z = 1;
  }

  // Initialize control point spacing
  _dx = (x2 - x1) / (_x - 1);
  _dy = (y2 - y1) / (_y - 1);
  if (z2 > z1) {
    _dz = (z2 - z1) / (_z - 1);
  } else {
    _dz = 1;
  }


  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Intialize memory for control point values
  _data = this->Allocate(_data, _x, _y, _z);

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z];
  if (_z == 1) {
    for (i = 0; i < 2*_x*_y*_z; i++) {
      _status[i] = _Active;
    }
    for (i = 2*_x*_y*_z; i < 3*_x*_y*_z; i++) {
      _status[i] = _Passive;
    }
  } else {
    for (i = 0; i < 3*_x*_y*_z; i++) {
      _status[i] = _Active;
    }
  }
}

irtkLinearFreeFormTransformation::irtkLinearFreeFormTransformation(irtkGenericImage<double> &image)
{
  int i, j, k;
  double x1, y1, z1, x2, y2, z2, xsize, ysize, zsize;

  // Check if image is displacement field
  if (image.GetT() != 3) {
    cerr << "irtkLinearFreeFormTransformation::irtkLinearFreeFormTransformation: Image does not contain displacement vectors";
    exit(1);
  }

  // Determine voxel size
  image.GetPixelSize(&xsize, &ysize, &zsize);

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
  _x = round((x2 - x1) / xsize) + 1;
  _y = round((y2 - y1) / ysize) + 1;
  if (z2 > z1) {
    _z = round((z2 - z1) / zsize) + 1;
  } else {
    _z = 1;
  }

  // Initialize control point spacing
  _dx = (x2 - x1) / (_x - 1);
  _dy = (y2 - y1) / (_y - 1);
  if (z2 > z1) {
    _dz = (z2 - z1) / (_z - 1);
  } else {
    _dz = 1;
  }

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Intialize memory for control point values
  _data = this->Allocate(_data, _x, _y, _z);

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z];
  if (_z == 1) {
    for (i = 0; i < 2*_x*_y*_z; i++) {
      _status[i] = _Active;
    }
    for (i = 2*_x*_y*_z; i < 3*_x*_y*_z; i++) {
      _status[i] = _Passive;
    }
  } else {
    for (i = 0; i < 3*_x*_y*_z; i++) {
      _status[i] = _Active;
    }
  }

  for (k = 0; k < image.GetZ(); k++) {
    for (j = 0; j < image.GetY(); j++) {
      for (i = 0; i < image.GetX(); i++) {
        this->Put(i, j, k, image(i, j, k, 0), image(i, j, k, 1), image(i, j, k, 2));
      }
    }
  }

}

irtkLinearFreeFormTransformation::irtkLinearFreeFormTransformation(double x1, double y1, double z1,
    double x2, double y2, double z2,
    double dx, double dy, double dz,
    double* xaxis, double* yaxis, double* zaxis)
{
  int i;

  // Initialize control point domain
  _origin._x = (x2 + x1) / 2.0;
  _origin._y = (y2 + y1) / 2.0;
  _origin._z = (z2 + z1) / 2.0;

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
  if (z2 > z1) {
    _z = round((z2 - z1) / dz) + 1;
  } else {
    _z = 1;
  }

  // Initialize control point spacing
  _dx = (x2 - x1) / (_x - 1);
  _dy = (y2 - y1) / (_y - 1);
  if (z2 > z1) {
    _dz = (z2 - z1) / (_z - 1);
  } else {
    _dz = 1;
  }

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Intialize memory for control point values
  _data = this->Allocate(_data, _x, _y, _z);

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z];
  if (_z == 1) {
    for (i = 0; i < 2*_x*_y*_z; i++) {
      _status[i] = _Active;
    }
    for (i = 2*_x*_y*_z; i < 3*_x*_y*_z; i++) {
      _status[i] = _Passive;
    }
  } else {
    for (i = 0; i < 3*_x*_y*_z; i++) {
      _status[i] = _Active;
    }
  }
}

irtkLinearFreeFormTransformation::irtkLinearFreeFormTransformation(const irtkLinearFreeFormTransformation &ffd) : irtkFreeFormTransformation3D(ffd)
{
  int i, j, k;

  // Initialize origin
  _origin = ffd._origin;

  // Initialize control point dimensions
  _x = ffd._x;
  _y = ffd._y;
  _z = ffd._z;

  // Intialize control point spacing
  _dx = ffd._dx;
  _dy = ffd._dy;
  _dz = ffd._dz;

  // Initialize x-axis
  _xaxis[0] = ffd._xaxis[0];
  _xaxis[1] = ffd._xaxis[1];
  _xaxis[2] = ffd._xaxis[2];

  // Initialize y-axis
  _yaxis[0] = ffd._yaxis[0];
  _yaxis[1] = ffd._yaxis[1];
  _yaxis[2] = ffd._yaxis[2];

  // Initialize y-axis
  _zaxis[0] = ffd._zaxis[0];
  _zaxis[1] = ffd._zaxis[1];
  _zaxis[2] = ffd._zaxis[2];

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Initialize memory for control point values
  _data = this->Allocate(_data, _x, _y, _z);
  for (k = -4; k < _z+4; k++) {
    for (j = -4; j < _y+4; j++) {
      for (i = -4; i < _x+4; i++) {
        _data[k][j][i]._x = ffd._data[k][j][i]._x;
        _data[k][j][i]._y = ffd._data[k][j][i]._y;
        _data[k][j][i]._z = ffd._data[k][j][i]._z;
      }
    }
  }

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z];
  for (i = 0; i < 3*_x*_y*_z; i++) {
    _status[i] = ffd._status[i];
  }
}

irtkLinearFreeFormTransformation::irtkLinearFreeFormTransformation(const irtkBSplineFreeFormTransformation &ffd)
{
  int i, j, k;
  double u, v, w, x, y, z;

  // Initialize origin
  _origin = ffd._origin;

  // Initialize control point dimensions
  _x = ffd._x;
  _y = ffd._y;
  _z = ffd._z;

  // Intialize control point spacing
  _dx = ffd._dx;
  _dy = ffd._dy;
  _dz = ffd._dz;

  // Initialize x-axis
  _xaxis[0] = ffd._xaxis[0];
  _xaxis[1] = ffd._xaxis[1];
  _xaxis[2] = ffd._xaxis[2];

  // Initialize y-axis
  _yaxis[0] = ffd._yaxis[0];
  _yaxis[1] = ffd._yaxis[1];
  _yaxis[2] = ffd._yaxis[2];

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Initialize memory for control point values
  _data = this->Allocate(_data, _x, _y, _z);
  for (k = -4; k < _z+4; k++) {
    for (j = -4; j < _y+4; j++) {
      for (i = -4; i < _x+4; i++) {
        x = i;
        y = j;
        z = k;
        u = x;
        v = y;
        w = z;
        // Calculate FFD
        ffd.FFD3D(u, v, w);

        // Calculate x, y, z
        ffd.LatticeToWorld(x, y, z);
        _data[k][j][i]._x = x + u;
        _data[k][j][i]._y = y + v;
        _data[k][j][i]._z = z + w;
      }
    }
  }

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z];
  for (i = 0; i < 3*_x*_y*_z; i++) {
    _status[i] = ffd._status[i];
  }
}

irtkLinearFreeFormTransformation::~irtkLinearFreeFormTransformation()
{
  // Free memory for control points if necessary
  if (_data != NULL) _data = this->Deallocate(_data, _x, _y, _z);

  _x = 0;
  _y = 0;
  _z = 0;
}

double irtkLinearFreeFormTransformation::Approximate(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no)
{
    int i, j, k, l, m, n, I, J, K, index;
    double s[2], t[2], u[2], x, y, z, B_I, B_J, B_K, basis;

    // Allocate memory for control points
    irtkVector3D<double> ***data = NULL;
    data = Allocate(data, _x, _y, _z);

    // Allocate memory for temporary storage
    double ***ds = NULL;
    ds = Allocate(ds, _x, _y, _z);

    // Initialize data structures
    for (k = -2; k < _z+2; k++) {
        for (j = -2; j < _y+2; j++) {
            for (i = -2; i < _x+2; i++) {
                data[k][j][i]._x = 0;
                data[k][j][i]._y = 0;
                data[k][j][i]._z = 0;
                ds[k][j][i] = 0;
            }
        }
    }
    // Initial loop: Calculate change of control points
    for (index = 0; index < no; index++) {
        x = x1[index];
        y = y1[index];
        z = z1[index];
        this->WorldToLattice(x, y, z);
        l = (int)floor(x);
        m = (int)floor(y);
        n = (int)floor(z);
        s[1] = x-l;
        t[1] = y-m;
        u[1] = z-n;
        s[0] = 1 - s[1];
        t[0] = 1 - t[1];
        u[0] = 1 - u[1];
        for (k = 0; k < 2; k++) {
            K = k + n;
            if ((K >= 0) && (K < _z)) {
                B_K = u[k];
                for (j = 0; j < 2; j++) {        
                    J = j + m;
                    if ((J >= 0) && (J < _y)) {
                        B_J = B_K*t[j];
                        for (i = 0; i < 2; i++) {
                            I = i + l;
                            if ((I >= 0) && (I < _x)) {
                                B_I = B_J * s[i];
                                basis = B_I;
                                data[K][J][I]._x += x2[index] * basis;
                                data[K][J][I]._y += y2[index] * basis;
                                data[K][J][I]._z += z2[index] * basis;
                                ds[K][J][I] += basis;
                            }
                        }
                    }
                }
            }
        }
    }

    // Final loop: Calculate new control points
    for (k = -2; k < _z+2; k++) {
        for (j = -2; j < _y+2; j++) {
            for (i = -2; i < _x+2; i++) {
                if (ds[k][j][i] > 0) {
                    _data[k][j][i]._x = data[k][j][i]._x / ds[k][j][i];
                    _data[k][j][i]._y = data[k][j][i]._y / ds[k][j][i];
                    _data[k][j][i]._z = data[k][j][i]._z / ds[k][j][i];
                }
            }
        }
    }

    // Calculate residual error
    for (index = 0; index < no; index++) {
        x = x1[index];
        y = y1[index];
        z = z1[index];
        this->LocalDisplacement(x, y, z);
        x2[index] -= x;
        y2[index] -= y;
        z2[index] -= z;
    }

    // Deallocate memory
    Deallocate(data, _x, _y, _z);
    Deallocate(ds, _x, _y, _z);

    return 0;
}

void irtkLinearFreeFormTransformation::ApproximateGradient(const double *x1, const double *y1, const double *z1, double *x2, double *y2, double *z2, int no, double *gradient)
{
    int i, j, k, l, m, n, I, J, K, index;
    double s[2], t[2], u[2], x, y, z, B_I, B_J, B_K, basis;

    // Allocate memory for control points
    irtkVector3D<double> ***data = NULL;
    data = Allocate(data, _x, _y, _z);

    // Allocate memory for temporary storage
    double ***ds = NULL;
    ds = Allocate(ds, _x, _y, _z);

    // Initialize data structures
    for (k = -2; k < _z+2; k++) {
        for (j = -2; j < _y+2; j++) {
            for (i = -2; i < _x+2; i++) {
                data[k][j][i]._x = 0;
                data[k][j][i]._y = 0;
                data[k][j][i]._z = 0;
                ds[k][j][i] = 0;
            }
        }
    }
    // Initial loop: Calculate change of control points
    for (index = 0; index < no; index++) {
        x = x1[index];
        y = y1[index];
        z = z1[index];
        this->WorldToLattice(x, y, z);
        l = (int)floor(x);
        m = (int)floor(y);
        n = (int)floor(z);
        s[1] = x-l;
        t[1] = y-m;
        u[1] = z-n;
        s[0] = 1 - s[1];
        t[0] = 1 - t[1];
        u[0] = 1 - u[1];
        for (k = 0; k < 2; k++) {
            K = k + n;
            if ((K >= 0) && (K < _z)) {
                B_K = u[k];
                for (j = 0; j < 2; j++) {        
                    J = j + m;
                    if ((J >= 0) && (J < _y)) {
                        B_J = B_K*t[j];
                        for (i = 0; i < 2; i++) {
                            I = i + l;
                            if ((I >= 0) && (I < _x)) {
                                B_I = B_J * s[i];
                                basis = B_I;
                                data[K][J][I]._x += x2[index] * basis;
                                data[K][J][I]._y += y2[index] * basis;
                                data[K][J][I]._z += z2[index] * basis;
                                ds[K][J][I] += basis;
                            }
                        }
                    }
                }
            }
        }
    }

    // Final loop: Calculate new control points
    for (k = 0; k < _z; k++) {
        for (j = 0; j < _y; j++) {
            for (i = 0; i < _x; i++) {
                index = this->LatticeToIndex(i,j,k);
                gradient[index] = data[k][j][i]._x;
                gradient[index + _x*_y*_z] = data[k][j][i]._y;
                gradient[index + 2*_x*_y*_z] = data[k][j][i]._z;
            }
        }
    }

    // Deallocate memory
    Deallocate(data, _x, _y, _z);
    Deallocate(ds, _x, _y, _z);
}

void irtkLinearFreeFormTransformation::Interpolate(const double *, const double *, const double *)
{
  cerr << "irtkLinearFreeFormTransformation::Interpolate: Not yet implemented" << endl;
  exit(1);
}

void irtkLinearFreeFormTransformation::Subdivide()
{
  cerr << "irtkLinearFreeFormTransformation::Subdivide: Not yet implemented" << endl;
  exit(1);
}

void irtkLinearFreeFormTransformation::FFD1(double &x, double &y, double &z) const
{
  int i, j, k;
  double t1, t2, u1, u2, v1, v2;

  // Check if there is some work to do
  if ((x < -1) || (y < -1) || (z < -1) || (x > _x) || (y > _y) || (z > _z)) {
    x = 0;
    y = 0;
    z = 0;
    return;
  }

  // Calculated integer coordinates
  i  = (int)floor(x);
  j  = (int)floor(y);
  k  = (int)floor(z);

  // Calculated fractional coordinates
  t1 = x - i;
  u1 = y - j;
  v1 = z - k;
  t2 = 1 - t1;
  u2 = 1 - u1;
  v2 = 1 - v1;

  // and use this...
  x = _data[k][j][i]._x;
  y = _data[k][j][i]._y;
  z = _data[k][j][i]._z;

  // Linear interpolation
  x = (t1 * (u2 * (v2 * _data[k][j][i+1]._x + v1 * _data[k+1][j][i+1]._x) +
             u1 * (v2 * _data[k][j+1][i+1]._x + v1 * _data[k+1][j+1][i+1]._x)) +
       t2 * (u2 * (v2 * _data[k][j][i]._x + v1 * _data[k+1][j][i]._x) +
             u1 * (v2 * _data[k][j+1][i]._x + v1 * _data[k+1][j+1][i]._x)));
  y = (t1 * (u2 * (v2 * _data[k][j][i+1]._y + v1 * _data[k+1][j][i+1]._y) +
             u1 * (v2 * _data[k][j+1][i+1]._y + v1 * _data[k+1][j+1][i+1]._y)) +
       t2 * (u2 * (v2 * _data[k][j][i]._y + v1 * _data[k+1][j][i]._y) +
             u1 * (v2 * _data[k][j+1][i]._y + v1 * _data[k+1][j+1][i]._y)));
  z = (t1 * (u2 * (v2 * _data[k][j][i+1]._z + v1 * _data[k+1][j][i+1]._z) +
             u1 * (v2 * _data[k][j+1][i+1]._z + v1 * _data[k+1][j+1][i+1]._z)) +
       t2 * (u2 * (v2 * _data[k][j][i]._z + v1 * _data[k+1][j][i]._z) +
             u1 * (v2 * _data[k][j+1][i]._z + v1 * _data[k+1][j+1][i]._z)));
}

void irtkLinearFreeFormTransformation::FFD2(double &x, double &y, double &z) const
{
  int i, j, k;
  double t1, t2, u1, u2, v1, v2;

  // Calculated integer coordinates
  i  = int(x);
  j  = int(y);
  k  = int(z);

  // Calculated fractional coordinates
  t1 = x - i;
  u1 = y - j;
  v1 = z - k;
  t2 = 1 - t1;
  u2 = 1 - u1;
  v2 = 1 - v1;

  // and use this...
  x = _data[k][j][i]._x;
  y = _data[k][j][i]._y;
  z = _data[k][j][i]._z;

  // Linear interpolation
  x = (t1 * (u2 * (v2 * _data[k][j][i+1]._x + v1 * _data[k+1][j][i+1]._x) +
             u1 * (v2 * _data[k][j+1][i+1]._x + v1 * _data[k+1][j+1][i+1]._x)) +
       t2 * (u2 * (v2 * _data[k][j][i]._x + v1 * _data[k+1][j][i]._x) +
             u1 * (v2 * _data[k][j+1][i]._x + v1 * _data[k+1][j+1][i]._x)));
  y = (t1 * (u2 * (v2 * _data[k][j][i+1]._y + v1 * _data[k+1][j][i+1]._y) +
             u1 * (v2 * _data[k][j+1][i+1]._y + v1 * _data[k+1][j+1][i+1]._y)) +
       t2 * (u2 * (v2 * _data[k][j][i]._y + v1 * _data[k+1][j][i]._y) +
             u1 * (v2 * _data[k][j+1][i]._y + v1 * _data[k+1][j+1][i]._y)));
  z = (t1 * (u2 * (v2 * _data[k][j][i+1]._z + v1 * _data[k+1][j][i+1]._z) +
             u1 * (v2 * _data[k][j+1][i+1]._z + v1 * _data[k+1][j+1][i+1]._z)) +
       t2 * (u2 * (v2 * _data[k][j][i]._z + v1 * _data[k+1][j][i]._z) +
             u1 * (v2 * _data[k][j+1][i]._z + v1 * _data[k+1][j+1][i]._z)));
}

void irtkLinearFreeFormTransformation::Jacobian(irtkMatrix &jac, double x, double y, double z, double t)
{
  this->LocalJacobian(jac, x, y, z, t);
}

void irtkLinearFreeFormTransformation::GlobalJacobian(irtkMatrix &jac, double, double, double, double)
{
  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Set matrix to identity
  jac(0, 0) = 1;
  jac(1, 1) = 1;
  jac(2, 2) = 1;
}

void irtkLinearFreeFormTransformation::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double)
{
  int i, j, k;
  double x1, y1, z1, x2, y2, z2;

  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);

  // Check if there is some work to do
  if ((x < 1) || (y < 1) || (z < 1) || (x >= _x-1) || (y >= _y-1) || (z >= _z-1)) {
    jac(0, 0) += 1;
    jac(1, 1) += 1;
    jac(2, 2) += 1;
    return;
  }

  // Compute derivatives
  i = (int)floor(x);
  j = (int)floor(y);
  k = (int)floor(z);

  // Jacobian matrix is 3 x 3
  jac.Initialize(3, 3);

  // Compute derivative in x
  x1 = x+0.5;
  y1 = y;
  z1 = z;
  x2 = x-0.5;
  y2 = y;
  z2 = z;
  this->FFD1(x1, y1, z1);
  this->FFD1(x2, y2, z2);
  jac(0, 0) = x1 - x2;
  jac(1, 0) = y1 - y2;
  jac(2, 0) = z1 - z2;

  // Compute derivative in y
  x1 = x;
  y1 = y+0.5;
  z1 = z;
  x2 = x;
  y2 = y-0.5;
  z2 = z;
  this->FFD1(x1, y1, z1);
  this->FFD1(x2, y2, z2);
  jac(0, 1) = x1 - x2;
  jac(1, 1) = y1 - y2;
  jac(2, 1) = z1 - z2;

  // Compute derivative in z
  x1 = x;
  y1 = y;
  z1 = z+0.5;
  x2 = x;
  y2 = y;
  z2 = z-0.5;
  this->FFD1(x1, y1, z1);
  this->FFD1(x2, y2, z2);
  jac(0, 2) = x1 - x2;
  jac(1, 2) = y1 - y2;
  jac(2, 2) = z1 - z2;

  jac.Transpose();
  // Convert derivatives to world coordinates
  jac = jac * _matW2L(0, 0, 3, 3);
  jac(0, 0) += 1;
  jac(1, 1) += 1;
  jac(2, 2) += 1;

}

double irtkLinearFreeFormTransformation::Bending(double, double, double, double
)
{
  cerr << "irtkLinearFreeFormTransformation::Bending :Not yet implemented" << endl;
  exit(1);
}

double irtkLinearFreeFormTransformation::Bending()
{
    int i, j, k;
    double bending;

    bending = 0;
    for (k = 0; k < _z; k++) {
        for (j = 0; j < _y; j++) {
            for (i = 0; i < _x; i++) {
                bending += this->Bending3D(i, j, k);
            }
        }
    }
    return bending;
}

void irtkLinearFreeFormTransformation::BendingGradient3D(double *gradient)
{
    // derivative of classic pairwise bending energy
    int index,i,j,k,x1,x2,y1,y2,z1,z2,n = _x*_y*_z;
    double tmp[3];
    for (k = 0; k < _z; k++) {
        for (j = 0; j < _y; j++) {
            for (i = 0; i < _x; i++) {
                index = this->LatticeToIndex(i, j, k);

                // evaluate gradient of the current point for bending energy

                x1 = i - 1;
                if(x1 < 0) x1 = 0;
                y1 = j - 1;
                if(y1 < 0) y1 = 0;
                z1 = k - 1;
                if(z1 < 0) z1 = 0;
                x2 = i + 1;
                if(x2 > this->_x - 1) x2 = this->_x - 1;
                y2 = j + 1;
                if(y2 > this->_y - 1) y2 = this->_y - 1;
                z2 = k + 1;
                if(z2 > this->_z - 1) z2 = this->_z - 1;

                tmp[0] = (2.0 * _data[k][j][i]._x - _data[z1][j][i]._x - _data[z2][j][i]._x)/_dz/_dz
                    + (2.0 * _data[k][j][i]._x - _data[k][y1][i]._x - _data[k][y2][i]._x)/_dy/_dy
                    + (2.0 * _data[k][j][i]._x - _data[k][j][x1]._x - _data[k][j][x2]._x)/_dx/_dx;
                tmp[1] = (2.0 * _data[k][j][i]._y - _data[z1][j][i]._y - _data[z2][j][i]._y)/_dz/_dz
                    + (2.0 * _data[k][j][i]._y - _data[k][y1][i]._y - _data[k][y2][i]._y)/_dy/_dy
                    + (2.0 * _data[k][j][i]._y - _data[k][j][x1]._y - _data[k][j][x2]._y)/_dx/_dx;
                tmp[2] = (2.0 * _data[k][j][i]._z - _data[z1][j][i]._z - _data[z2][j][i]._z)/_dz/_dz
                    + (2.0 * _data[k][j][i]._z - _data[k][y1][i]._z - _data[k][y2][i]._z)/_dy/_dy
                    + (2.0 * _data[k][j][i]._z - _data[k][j][x1]._z - _data[k][j][x2]._z)/_dx/_dx;

                gradient[index]     = -2.0*tmp[0];
                gradient[index+n]   = -2.0*tmp[1];
                gradient[index+2*n] = -2.0*tmp[2];

                // evaluate tested by numerical! good

            }
        }
    }
}

void irtkLinearFreeFormTransformation::BendingGradient(double *gradient)
{
    this->BendingGradient3D(gradient);
}

double irtkLinearFreeFormTransformation::Bending3D(int i, int j, int k)
{
    // classic pairwise bending energy
    double dx, dy, dz;
    int x2,y2,z2;

    x2 = i + 1;
    if(x2 > this->_x - 1) x2 = this->_x - 1;
    y2 = j + 1;
    if(y2 > this->_y - 1) y2 = this->_y - 1;
    z2 = k + 1;
    if(z2 > this->_z - 1) z2 = this->_z - 1;

    if (i != x2) {
        dx = pow((_data[k][j][x2]._x - _data[k][j][i]._x)/(x2 - i)/_dx,2.0)
            + pow((_data[k][j][x2]._y - _data[k][j][i]._y)/(x2 - i)/_dx,2.0)
            + pow((_data[k][j][x2]._z - _data[k][j][i]._z)/(x2 - i)/_dx,2.0);
    }else{
        dx = 0;
    }

    if (j != y2) {
        dy = pow((_data[k][y2][i]._x - _data[k][j][i]._x)/(y2 - j)/_dy,2.0)
            + pow((_data[k][y2][i]._y - _data[k][j][i]._y)/(y2 - j)/_dy,2.0)
            + pow((_data[k][y2][i]._z - _data[k][j][i]._z)/(y2 - j)/_dy,2.0);
    }else{
        dy = 0;
    }

    if (k != z2) {
        dz = pow((_data[z2][j][i]._x - _data[k][j][i]._x)/(z2 - k)/_dz,2.0)
            + pow((_data[z2][j][i]._y - _data[k][j][i]._y)/(z2 - k)/_dz,2.0)
            + pow((_data[z2][j][i]._z - _data[k][j][i]._z)/(z2 - k)/_dz,2.0);
    }else{
        dz = 0;
    }

    return dx + dy + dz;
}

void irtkLinearFreeFormTransformation::BoundingBoxCP(int index, irtkPoint &p1, irtkPoint &p2, double fraction) const
{
  int i, j, k;

  if (index >= _x*_y*_z) {
    index -= _x*_y*_z;
    if (index >= _x*_y*_z) {
      index -= _x*_y*_z;
    }
  }
  i = index/(_y*_z);
  j = index%(_y*_z)/_z;
  k = index%(_y*_z)%_z;
  p1 = irtkPoint(i-fraction, j-fraction, k-fraction);
  p2 = irtkPoint(i+fraction, j+fraction, k+fraction);
  this->LatticeToWorld(p1);
  this->LatticeToWorld(p2);
}

void irtkLinearFreeFormTransformation::BoundingBoxImage(irtkGreyImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
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

void irtkLinearFreeFormTransformation::BoundingBoxImage(irtkRealImage *image, int index, int &i1, int &j1, int &k1, int &i2, int &j2, int &k2, double fraction) const
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

int irtkLinearFreeFormTransformation::CheckHeader(char *name)
{
  char buffer[255];

  ifstream from(name);

  if (!from) {
    cerr << "irtkLinearFreeFormTransformation::CheckHeader: Can't open file "
    << name << "\n";
    exit(1);
  }

  // Read keyword
  from >> buffer;
  if (strcmp(buffer, "LinearFFD:") != 0) {
    return false;
  }

  return true;
}

istream& irtkLinearFreeFormTransformation::Import(istream& is)
{
  char buffer[255];
  int i, j, k, n, index, no_attr;
  double x1, y1, z1, x2, y2, z2;

  // Free memory if necessary
  _data = Deallocate(_data, _x, _y, _z);

  delete []_status;
  _x = 0;
  _y = 0;
  _z = 0;

  // Read keyword
  is >> buffer;

  if (strcmp(buffer, "LinearFFD:") == 0) {

    // Read no. of control points
    is >> buffer;
    is >> buffer;
    is >> buffer >> _x  >> buffer >> _y
    >> buffer >> _z;

    // Read bounding box for control points
    is >> buffer >> x1 >> y1 >> z1 >> buffer >> x2 >> y2 >> z2;

    // Calculate origin
    _origin._x = (x2 + x1) / 2.0;
    _origin._y = (y2 + y1) / 2.0;
    _origin._z = (z2 + z1) / 2.0;

    // Calculate control point spacing
    _dx = (x2 - x1) / (_x-1) ;
    _dy = (y2 - y1) / (_y-1) ;
    if (z2 > z1) {
      _dz = (z2 - z1) / (_z-1);
    } else {
      _dz = 1;
    }

    // Skip rest of line
    is.get(buffer, 255);
    is.clear();
    is.seekg(1, ios::cur);

  } else {
    if (strcmp(buffer, "LinearFFD_BINARY:") == 0) {

      // Read no. of control points
      is >> buffer;
      is >> buffer;
      is >> buffer >> _x  >> buffer >> _y >> buffer >> _z;

      // Skip rest of line
      is.get(buffer, 255);
      is.clear();
      is.seekg(1, ios::cur);

      // Read orientation, origin and spacing
      double tmp;

      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _xaxis[0] = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _xaxis[1] = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _xaxis[2] = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _yaxis[0] = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _yaxis[1] = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _yaxis[2] = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _dx = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _dy = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _dz = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _origin._x = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _origin._y = tmp;
      is.read((char *)&tmp, sizeof(double));
#ifndef WORDS_BIGENDIAN
      swap64((char *)&tmp, (char *)&tmp, 1);
#endif
      _origin._z = tmp;

    } else {
      cerr << "istream& operator>> (istream& is, irtkLinearFreeFormTransformation &transformation): Invalid file format: " << buffer << endl;
      exit(1);
    }
  }

  // Compute z-axis.
  _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
  _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
  _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];

  // Calculate number of control points
  n = _x*_y*_z;

  // Initialize control points
  _data = Allocate(_data, _x, _y, _z);

  // Initialize control points
  for (i = -2; i < _x+2; i++) {
    for (j = -2; j < _y+2; j++) {
      for (k = -2; k < _z+2; k++) {
        _data[k][j][i]._x = 0;
        _data[k][j][i]._y = 0;
        _data[k][j][i]._z = 0;
      }
    }
  }

  // Initialize memory for control point status
  _status = new _Status[3*n];

  // Initialize control point status to active
  for (i = 0; i < 3*n; i++) {
    _status[i] = _Active;
  }

  // Allocate temporary memory
  float *data = new float[3*n*sizeof(float)];

  // Read binary data for control point values
  is.read((char *)data, 3*n*sizeof(float));

#ifndef WORDS_BIGENDIAN
  swap32((char *)data, (char *)data, 3*n);
#endif

  // Convert data
  index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        _data[k][j][i]._x = data[index];
        _data[k][j][i]._y = data[index+1];
        _data[k][j][i]._z = data[index+2];
        index += 3;
      }
    }
  }

  // Free temporary memory
  delete []data;

  // Skip rest of line
  is.get(buffer, 255);
  is.clear();
  is.seekg(1, ios::cur);

  // Read number of attributes
  is >> buffer;

  if (strcmp(buffer, "ATTRIBUTES:") == 0) {
    is >> no_attr;
  } else {
    cerr << "istream& operator>> irtkLinearFreeFormTransformation: "
    << "Expecting attributes " << buffer << endl;
    exit(1);
  }

  for (i = 0; i < no_attr; i++) {

    // Read attribute
    is >> buffer;

    // Allocate temporary memory
    unsigned char *status = new unsigned char[n*sizeof(unsigned char)];

    if (strcmp(buffer, "STATUS_X:") == 0) {

      // Skip rest of line
      is.get(buffer, 255);
      is.clear();
      is.seekg(1, ios::cur);
      // Read binary data for control point status
      is.read((char *)&(status[0]), n*sizeof(unsigned char));

      // Convert status
      for (j = 0; j < n; j++) {
        _status[j] = (_Status) status[j];
      }
    } else {
      if (strcmp(buffer, "STATUS_Y:") == 0) {

        // Skip rest of line
        is.get(buffer, 255);
        is.clear();
        is.seekg(1, ios::cur);
        // Read binary data for control point status
        is.read((char *)&(status[0]), n*sizeof(unsigned char));

        // Convert status
        for (j = 0; j < n; j++) {
          _status[j+n] = (_Status) status[j];
        }
      } else {
        if (strcmp(buffer, "STATUS_Z:") == 0) {

          // Skip rest of line
          is.get(buffer, 255);
          is.clear();
          is.seekg(1, ios::cur);
          // Read binary data for control point status
          is.read((char *)&(status[0]), n*sizeof(unsigned char));

          // Convert status
          for (j = 0; j < n; j++) {
            _status[j+2*n] = (_Status) status[j];
          }
        } else {
          cerr << "istream& operator>> irtkLinearFreeFormTransformation: ";
          cerr << "Unknown attribute:" << buffer << " terminating" << endl;
          exit(1);
        }
      }
    }

    // Free temporary memory
    delete [] status;

    // Skip rest of line
    is.get(buffer, 255);
    is.clear();
    is.seekg(1, ios::cur);

  }
  // Update transformation matrix
  this->UpdateMatrix();

  return is;
}

ostream& irtkLinearFreeFormTransformation::Export(ostream& os)
{
  int i, j, k, n, index;
  double x1, y1, z1, x2, y2, z2;
  double zaxis[3];

  cerr << "irtkLinearFreeFormTransformation::Export: DO NOT USE. THIS METHOD IS OBSOLETE" << endl;

  // Compute z-axis.
  zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
  zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
  zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];

  if ((zaxis [0] != _zaxis[0]) || (zaxis [1] != _zaxis[1]) || (zaxis [2] != _zaxis[2])) {
    cerr << "irtkLinearFreeFormTransformation::Export: z-axis has not default value. Writing transformation will omit z-axis." << endl;
  }

  // Check if we need to save orientation
  if ((_xaxis[0] != 1) || (_xaxis[1] != 0) ||
      (_xaxis[2] != 0) || (_yaxis[0] != 0) ||
      (_yaxis[1] != 1) || (_yaxis[2] != 0)) {

    // Write keyword and no. of DOFs
    os << "LinearFFD_BINARY: " << this->NumberOfDOFs() << endl;

    // Write no. of control points
    os << "Control points: " << _x << " x " <<
    _y << " x " << _z << endl;

    // Write orientation, origin and spacing
    double tmp;

    tmp = _xaxis[0];
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _xaxis[1];
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _xaxis[2];
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _yaxis[0];
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _yaxis[1];
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _yaxis[2];
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _dx;
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _dy;
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _dz;
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _origin._x;
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _origin._y;
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));
    tmp = _origin._z;
#ifndef WORDS_BIGENDIAN
    swap64((char *)&tmp, (char *)&tmp, 1);
#endif
    os.write((char *)&tmp, sizeof(double));

  } else {

    // Write keyword and no. of DOFs
    os << "LinearFFD: " << this->NumberOfDOFs() << endl;

    // Write no. of control points
    os << "Control points: " << _x << " x " <<
    _y << " x " << _z << endl;

    // Calculate bounding box
    x1 = 0;
    y1 = 0;
    z1 = 0;
    this->LatticeToWorld(x1, y1, z1);
    x2 = _x-1;
    y2 = _y-1;
    z2 = _z-1;
    this->LatticeToWorld(x2, y2, z2);

    // Write bounding box
    os << "Point_1: " << x1 << " " << y1 << " " << z1 << " ";
    os << "Point_2: " << x2 << " " << y2 << " " << z2 << endl;
  }

  // Compute no. of control points
  n = _x*_y*_z;

  // Allocate temporary memory
  float *data = new float[3*n*sizeof(float)];

  // Convert data
  index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        data[index]   = _data[k][j][i]._x;
        data[index+1] = _data[k][j][i]._y;
        data[index+2] = _data[k][j][i]._z;
        index += 3;
      }
    }
  }

#ifndef WORDS_BIGENDIAN
  swap32((char *)data, (char *)data, 3*n);
#endif

  // Write binary data for control point values
  os.write((char *)data, 3*n*sizeof(float));

  // Free temporary memory
  delete []data;

  // Write end of line
  os << endl;

  // Cast status (later label) to unsigned char to avoid byte order problems
  unsigned char *status = new unsigned char [3*n*sizeof(unsigned char)];
  for (index = 0; index < 3*n; index++) {
    status[index] = _status[index];
  }

  // Write attributes
  os << "ATTRIBUTES: 3" << endl;

  // Write data for control point status
  os << "STATUS_X:" << endl;
  os.write((char *)&(status[0]), n*sizeof(unsigned char));
  os << endl;
  os << "STATUS_Y:" << endl;
  os.write((char *)&(status[n]), n*sizeof(unsigned char));
  os << endl;
  os << "STATUS_Z:" << endl;
  os.write((char *)&(status[2*n]), n*sizeof(unsigned char));

  // Free temporary memory
  delete []status;

  // Write end of line
  os << endl;

  return os;
}

void irtkLinearFreeFormTransformation::Print()
{
  // Write keyword and no. of DOFs
  cout << "LinearFFD: " << this->NumberOfDOFs() << endl;

  // Write no. of control points
  cout << "Control points: " << _x << " x " << _y << " x " << _z << endl;
  cout << "Spacing: " << _dx << " x " << _dy << " x " << _dz << endl;
  cout << "Origin: " << _origin._x << " " << _origin._y << " " << _origin._z << " ";
  cout << "Orientation: " << _xaxis[0] << " " << _xaxis[1] << " "
  << _xaxis[2] << " " << _yaxis[0] << " " << _yaxis[1] << " "
  << _yaxis[2] << endl;
}

irtkCifstream& irtkLinearFreeFormTransformation::Read(irtkCifstream& from)
{
  double *data;
  int i, j, k, index;
  unsigned int magic_no, trans_type;

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1);
  if (magic_no != IRTKTRANSFORMATION_MAGIC) {
    cerr << "irtkLinearFreeFormTransformation::Read: Not a vaild transformation file" << endl;
    exit(1);
  }

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);
  if ((trans_type != IRTKTRANSFORMATION_LINEAR_FFD) && (trans_type != IRTKTRANSFORMATION_LINEAR_FFD_EXT1)) {
    cerr << "irtkLinearFreeFormTransformation::Read: Not a vaild linear FFD transformation" << endl;
    exit(1);
  }

  // Free memory if necessary
  _data = Deallocate(_data, _x, _y, _z);
  delete []_status;

  // Read no of control points
  from.ReadAsInt(&_x, 1);
  from.ReadAsInt(&_y, 1);
  from.ReadAsInt(&_z, 1);

  // Read orientation of bounding box
  from.ReadAsDouble(_xaxis, 3);
  from.ReadAsDouble(_yaxis, 3);
  if (trans_type == IRTKTRANSFORMATION_LINEAR_FFD_EXT1) {
    from.ReadAsDouble(_zaxis, 3);
  } else {
    _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
    _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
    _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];
  }

  // Read spacing of bounding box
  from.ReadAsDouble(&_dx, 1);
  from.ReadAsDouble(&_dy, 1);
  from.ReadAsDouble(&_dz, 1);

  // Read spacing of bounding box
  from.ReadAsDouble(&_origin._x, 1);
  from.ReadAsDouble(&_origin._y, 1);
  from.ReadAsDouble(&_origin._z, 1);

  // Initialize control points
  _data = Allocate(_data, _x, _y, _z);

  // Initialize control points
  for (i = -2; i < _x+2; i++) {
    for (j = -2; j < _y+2; j++) {
      for (k = -2; k < _z+2; k++) {
        _data[k][j][i]._x = 0;
        _data[k][j][i]._y = 0;
        _data[k][j][i]._z = 0;
      }
    }
  }

  // Allocate temporary memory
  data = new double[3*_x*_y*_z];

  // Read control point data
  from.ReadAsDouble(data, 3*_x*_y*_z);

  // Convert data
  index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        _data[k][j][i]._x = data[index];
        _data[k][j][i]._y = data[index+1];
        _data[k][j][i]._z = data[index+2];
        index += 3;
      }
    }
  }

  // Free temporary memory
  delete []data;

  // Initialize memory for control point status
  _status = new _Status[3*_x*_y*_z];

  // Read control point status
  from.ReadAsInt((int *)_status, 3*_x*_y*_z);

  // Update transformation matrix
  this->UpdateMatrix();

  return from;
}

irtkCofstream& irtkLinearFreeFormTransformation::Write(irtkCofstream& to)
{
  double *data;
  int i, j, k, index;
  unsigned int magic_no, trans_type;

  // Write magic no. for transformations
  magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type
  trans_type = IRTKTRANSFORMATION_LINEAR_FFD_EXT1;
  to.WriteAsUInt(&trans_type, 1);

  // Write no of control points
  to.WriteAsInt(&_x, 1);
  to.WriteAsInt(&_y, 1);
  to.WriteAsInt(&_z, 1);

  // Write orientation of bounding box
  to.WriteAsDouble(_xaxis, 3);
  to.WriteAsDouble(_yaxis, 3);
  to.WriteAsDouble(_zaxis, 3);

  // Write spacing of bounding box
  to.WriteAsDouble(&_dx, 1);
  to.WriteAsDouble(&_dy, 1);
  to.WriteAsDouble(&_dz, 1);

  // Write spacing of bounding box
  to.WriteAsDouble(&_origin._x, 1);
  to.WriteAsDouble(&_origin._y, 1);
  to.WriteAsDouble(&_origin._z, 1);

  // Allocate temporary memory
  data = new double[3*_x*_y*_z];

  // Convert data
  index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        data[index]   = _data[k][j][i]._x;
        data[index+1] = _data[k][j][i]._y;
        data[index+2] = _data[k][j][i]._z;
        index += 3;
      }
    }
  }

  // Write control point data
  to.WriteAsDouble(data, 3*_x*_y*_z);

  // Free temporary memory
  delete []data;

  // Write control point status
  to.WriteAsInt((int *)_status, 3*_x*_y*_z);

  return to;
}

double irtkLinearFreeFormTransformation::Inverse(double &x, double &y, double &z, double, double)
{
  /*
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
    float invert[3], f_invert[3];
    invert[0] = x;
    invert[1] = y;
    invert[2] = z;

    // Numerically approximate the inverse transformation
    newt2(invert-1, 3, &check, Newton_function);

    // Calculate error
    irtkTransformationEvaluate(3, invert-1, f_invert-1);
    error = sqrt(f_invert[0]*f_invert[0]+f_invert[1]*f_invert[1]+f_invert[2]*f_invert[2]);
    if (error > tolerance) {
      cout << "irtkFreeFormTransformation3D::Inverse: RMS error = " << error << "\n";
    }

    // Set output to solution
    x = invert[0];
    y = invert[1];
    z = invert[2];

    return error;
  */

  double u, v, w;

  u = x;
  v = y;
  w = z;

  // Convert world coordinates in to FFD coordinates
  this->WorldToLattice(u, v, w);

  // Calculate FFD
  this->FFD1(u, v, w);

  // Add FFD to world coordinates
  x -= u;
  y -= v;
  z -= w;

  return 0;
}

void irtkLinearFreeFormTransformation::Compose(irtkTransformation *t1)
{
  int i, j, k;
  double x1, y1, z1, x2, y2, z2;

  // Allocate transformation
  irtkLinearFreeFormTransformation *t2 = new irtkLinearFreeFormTransformation(*this);

  // Compose t2 o t1
  for (k = 0; k < this->GetZ(); k++) {
    for (j = 0; j < this->GetY(); j++) {
      for (i = 0; i < this->GetX(); i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        this->LatticeToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        // Apply first transformation
        t1->Transform(x2, y2, z2);
        // Apply second transformation
        t2->Transform(x2, y2, z2);
        // Update to displacement
        this->Put(i, j, k, x2 - x1, y2 - y1, z2 - z1);
      }
    }
  }

  // Delete transformation
  delete t2;
}
