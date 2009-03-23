/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkBiasField.h>

#include <nr.h>
#include <nrutil.h>

#define LUTSIZE (double)(BIASLOOKUPTABLESIZE-1)

double irtkBSplineBiasField::LookupTable   [BIASLOOKUPTABLESIZE][4];

double irtkBSplineBiasField::LookupTable_I [BIASLOOKUPTABLESIZE][4];

double irtkBSplineBiasField::LookupTable_II[BIASLOOKUPTABLESIZE][4];

irtkBSplineBiasField::irtkBSplineBiasField()
{
  int i;

  // Initialize control point domain
  _origin._x = 0;
  _origin._y = 0;
  _origin._z = 0;

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

irtkBSplineBiasField::irtkBSplineBiasField(const irtkGreyImage &image, double dx, double dy, double dz)
{
  int i;
  double x1, y1, z1, x2, y2, z2;

  cerr<<"Initialising bias field using control points spacing: "<< dx << " " << dy << " " << dz <<endl;

  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = image.GetX()-1;
  y2 = image.GetY()-1;
  z2 = image.GetZ()-1;
  image.ImageToWorld(x1, y1, z1);
  image.ImageToWorld(x2, y2, z2);
  image.GetOrientation(_xaxis, _yaxis, _zaxis);

  // Initialize control point domain
  _origin._x = (x2 + x1) / 2.0;
  _origin._y = (y2 + y1) / 2.0;
  _origin._z = (z2 + z1) / 2.0;

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
  if (x2 > x1) {
    _x = round((x2 - x1) / dx) + 1;
  } else {
    _x = 1;
  }
  if (y2 > y1) {
    _y = round((y2 - y1) / dy) + 1;
  } else {
    _y = 1;
  }
  if (z2 > z1) {
    _z = round((z2 - z1) / dz) + 1;
  } else {
    _z = 1;
  }

  // Initialize control point spacing
  if (x2 > x1) {
    _dx = (x2 - x1) / (_x - 1);
  } else {
    _dx = 1;
  }
  if (y2 > y1) {
    _dy = (y2 - y1) / (_y - 1);
  } else {
    _dy = 1;
  }
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

  // Initialize lookup table
  for (i = 0; i < BIASLOOKUPTABLESIZE; i++) {
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

irtkBSplineBiasField::irtkBSplineBiasField(const irtkGreyImage &image, int x, int y, int z, bool bounding_box, int padding )
{
//arguments: number of control points
  int i,j,k;
  double x1, y1, z1, x2, y2, z2;

  cerr<<"Initialising bias field using number of control points: "<< x << "x" << y << "x" << z <<endl;


  x1 = 0;
  y1 = 0;
  z1 = 0;
  x2 = image.GetX()-1;
  y2 = image.GetY()-1;
  z2 = image.GetZ()-1;

  if (bounding_box) {

    cerr<<"Bounding box: ";
    ///x1
    i=0; j=0; k=0;
    bool found = false;
    while ( (!found) && (i<image.GetX()) ) {
      while ( (!found) && (j<image.GetY())) {
        while ( (!found) && (k<image.GetZ())) {
          if (image.Get(i,j,k) > padding) found = true;
          k++;
        }
        j++;
      }
      i++;
    }
    x1 = i;
    cerr<<"found=" <<found <<", ";
    cerr<<"padding=" <<padding <<", ";
    cerr<<"x1=" <<x1 <<", ";

    ///x2
    i=image.GetX()-1; j=0; k=0;
    found = false;
    while ( (!found) && (i>0) ) {
      while ( (!found) && (j<image.GetY())) {
        while ( (!found) && (k<image.GetZ())) {
          if (image.Get(i,j,k) > padding) found = true;
          k++;
        }
        j++;
      }
      i--;
    }
    x2 = i;
    cerr<<"x2=" <<x2 <<", ";

  }

  image.ImageToWorld(x1, y1, z1);
  image.ImageToWorld(x2, y2, z2);
  image.GetOrientation(_xaxis, _yaxis, _zaxis);

  // Initialize control point domain
  _origin._x = (x2 + x1) / 2.0;
  _origin._y = (y2 + y1) / 2.0;
  _origin._z = (z2 + z1) / 2.0;

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
  _x=x;
  _y=y;
  _z=z;

  /*if (x2 > x1){
    _x = round((x2 - x1) / dx) + 1;
  } else {
    _x = 1;
  }
  if (y2 > y1){
    _y = round((y2 - y1) / dy) + 1;
  } else {
    _y = 1;
  }
  if (z2 > z1){
    _z = round((z2 - z1) / dz) + 1;
  } else {
    _z = 1;
  }
  */
  // Initialize control point spacing
  if (x2 > x1) {
    _dx = (x2 - x1) / (_x - 1);
  } else {
    _dx = 1;
  }
  if (y2 > y1) {
    _dy = (y2 - y1) / (_y - 1);
  } else {
    _dy = 1;
  }
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

  // Initialize lookup table
  for (i = 0; i < BIASLOOKUPTABLESIZE; i++) {
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

irtkBSplineBiasField::irtkBSplineBiasField(const irtkBSplineBiasField &bias)
{
  int i, j, k;

  // Initialize origin
  _origin = bias._origin;

  // Initialize control point dimensions
  _x = bias._x;
  _y = bias._y;
  _z = bias._z;

  // Intialize control point spacing
  _dx = bias._dx;
  _dy = bias._dy;
  _dz = bias._dz;

  // Initialize x-axis
  _xaxis[0] = bias._xaxis[0];
  _xaxis[1] = bias._xaxis[1];
  _xaxis[2] = bias._xaxis[2];

  // Initialize y-axis
  _yaxis[0] = bias._yaxis[0];
  _yaxis[1] = bias._yaxis[1];
  _yaxis[2] = bias._yaxis[2];

  // Initialize z-axis
  _zaxis[0] = bias._zaxis[0];
  _zaxis[1] = bias._zaxis[1];
  _zaxis[2] = bias._zaxis[2];

  // Initialize transformation matrix
  _matL2W = irtkMatrix(4, 4);
  _matW2L = irtkMatrix(4, 4);

  // Update transformation matrix
  this->UpdateMatrix();

  // Initialize memory for control point values
  _data = this->Allocate(_data, _x, _y, _z);
  for (i = -4; i < _x+4; i++) {
    for (j = -4; j < _y+4; j++) {
      for (k = -4; k < _z+4; k++) {
        _data[k][j][i] = bias._data[k][j][i];
      }
    }
  }

  // Initialize lookup table
  for (i = 0; i < BIASLOOKUPTABLESIZE; i++) {
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

irtkBSplineBiasField::~irtkBSplineBiasField()
{
  // Free memory for control points if necessary
  if (_data != NULL) _data = this->Deallocate(_data, _x, _y, _z);

  _x = 0;
  _y = 0;
  _z = 0;
}

double irtkBSplineBiasField::FFD1(double x, double y, double z) const
{
  // Check if there is some work to do
  if ((x < -2) || (y < -2) || (z < -2) || (x > _x+1) ||
      (y > _y+1) || (z > _z+1)) {
    return 0;
  }

  double *data;
  double s, t, u, B_I, B_J, B_K, xi, xii;
  int i, j, k, l, m, n, S, T, U;

  // Now calculate the real stuff
  l = (int)floor(x);
  m = (int)floor(y);
  n = (int)floor(z);
  s = x-l;
  t = y-m;
  u = z-n;

  // Calculate offset
  i = (_x + 8) * (_y + 4);
  x = 0;
  S = round(LUTSIZE*s);
  T = round(LUTSIZE*t);
  U = round(LUTSIZE*u);
  data = &(_data[n-1][m-1][l-1]);
  for (k = 0; k < 4; k++) {
    B_K = this->LookupTable[U][k];
    xi = 0;
    for (j = 0; j < 4; j++) {
      B_J = this->LookupTable[T][j];

      // Inner most loop unrolled starts here
      B_I = this->LookupTable[S][0];
      xii = *data * this->LookupTable[S][0];
      data++;
      B_I = this->LookupTable[S][1];
      xii += *data * B_I;
      data++;
      B_I = this->LookupTable[S][2];
      xii += *data * B_I;
      data++;
      B_I = this->LookupTable[S][3];
      xii += *data * B_I;
      data++;
      // Inner most loop unrolled stops here

      xi += xii * B_J;
      data += _x + 4;
    }
    x += xi * B_K;
    data += i;
  }
  return 1+x;
}

double irtkBSplineBiasField::FFD2(double x, double y, double z) const
{
  double *data;
  double s, t, u, B_I, B_J, B_K, xi, xii;
  int i, j, k, l, m, n, S, T, U;

  // Now calculate the real stuff
  l = (int)floor(x);
  m = (int)floor(y);
  n = (int)floor(z);
  s = x-l;
  t = y-m;
  u = z-n;

  // Calculate offset
  i = (_x + 8) * (_y + 4);
  x = 0;
  S = round(LUTSIZE*s);
  T = round(LUTSIZE*t);
  U = round(LUTSIZE*u);
  data = &(_data[n-1][m-1][l-1]);
  for (k = 0; k < 4; k++) {
    B_K = this->LookupTable[U][k];
    xi = 0;
    for (j = 0; j < 4; j++) {
      B_J = this->LookupTable[T][j];

      // Inner most loop unrolled starts here
      B_I = this->LookupTable[S][0];
      xii = *data * this->LookupTable[S][0];
      data++;
      B_I = this->LookupTable[S][1];
      xii += *data * B_I;
      data++;
      B_I = this->LookupTable[S][2];
      xii += *data * B_I;
      data++;
      B_I = this->LookupTable[S][3];
      xii += *data * B_I;
      data++;
      // Inner most loop unrolled stops here

      xi += xii * B_J;
      data += _x + 4;
    }
    x += xi * B_K;
    data += i;
  }
  return 1+x;
}

double irtkBSplineBiasField::Approximate(double *x1, double *y1, double *z1, double *bias, int no)
{
  int i, j, k, l, m, n, I, J, K, index;
  double b, s, t, u, x, y, z, error, phi, norm, tmp;

  // Allocate memory
  double ***dx = NULL;
  double ***ds = NULL;
  dx = ::Allocate(dx, _x, _y, _z);
  ds = ::Allocate(ds, _x, _y, _z);

  // Subtract displacements which are approximated by current control points
  for (index = 0; index < no; index++) {
    x = x1[index];
    y = y1[index];
    z = z1[index];
    bias[index] -= this->Bias(x, y, z);
  }

  // Initialize data structures
  for (k = 0; k < _z; k++) {
    for (j = 0; j < _y; j++) {
      for (i = 0; i < _x; i++) {
        dx[k][j][i] = 0;
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
    s = x-l;
    t = y-m;
    u = z-n;
    norm = 0;
    for (k = 0; k < 4; k++) {
      for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++) {
          tmp   = B(i, s) * B(j, t) * B(k, u);
          norm += tmp * tmp;
        }
      }
    }
    for (k = 0; k < 4; k++) {
      K = k + n - 1;
      if ((K >= 0) && (K < _z)) {
        for (j = 0; j < 4; j++) {
          J = j + m - 1;
          if ((J >= 0) && (J < _y)) {
            for (i = 0; i < 4; i++) {
              I = i + l - 1;
              if ((I >= 0) && (I < _x)) {
                b = B(i, s) * B(j, t) * B(k, u);
                phi = bias[index] * b / norm;
                dx[K][J][I] += b * b * phi;
                ds[K][J][I] += b * b;
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
    bias[index] += this->Bias(x, y, z);
  }

  // Final loop: Calculate new control points
  for (k = 0; k < _z; k++) {
    for (j = 0; j < _y; j++) {
      for (i = 0; i < _x; i++) {
        if (ds[k][j][i] > 0) {
          _data[k][j][i] += dx[k][j][i] / ds[k][j][i];
        }
      }
    }
  }

  // Calculate residual error
  error = 0;
  for (index = 0; index < no; index++) {
    x = x1[index];
    y = y1[index];
    z = z1[index];
    bias[index] -= this->Bias(x, y, z);

    // Calculate error
    error += sqrt(bias[index]*bias[index]);
  }
  error = error / (double)no;

  // Deallocate memory
  ::Deallocate<double>(dx);
  ::Deallocate<double>(ds);

  // Return error
  return error;
}

void irtkBSplineBiasField::WeightedLeastSquares(double *x1, double *y1, double *z1, double *bias, double *weights, int no)
{
  int index, i, j, k, ii, jj, kk, l, m, n, indx, indx1;
  double b,w, s, t, u, x, y, z;
  irtkVector P(_x*_y*_z),T(_x*_y*_z);
  irtkMatrix M(_x*_y*_z, _x*_y*_z);

  cerr<<_x<<" "<<_y<<" "<< _z<<endl;

  cerr<<"Starting weighted least squares...";

  int per=0;

  for (index = 0; index < no; index++) {
    if (index*10/no > per) {
      per++;
      cerr<<per<<"0%...";
    }
    //cerr<<index<<" "<<endl;
    x = x1[index];
    y = y1[index];
    z = z1[index];
    b = bias[index];
    w = weights[index];

    //cerr<<x<<" "<<y<<" "<<z<<endl;

    this->WorldToLattice(x, y, z);

    ///this must be added because of small numerical errors introduced by change between the coordinate systems
    if (x<0) x=0;
    if (y<0) y=0;
    if (z<0) z=0;
    if (x>_x-1) x=_x-1;
    if (y>_y-1) y=_y-1;
    if (z>_z-1) z=_z-1;

    l = (int)floor(x);
    m = (int)floor(y);
    n = (int)floor(z);
    s = x-l;
    t = y-m;
    u = z-n;

    ///Adjust parameters for x=_x,...
    if (x==_x-1) {
      l=_x-2;
      s=1;
    }
    if (y==_y-1) {
      m=_y-2;
      t=1;
    }

    if (z==_z-1) {
      n=_z-2;
      u=1;
    }

    //cerr<<s<<" "<<t<<" "<<u<<" "<<w<<" "<<b<<endl;


    for (k = 0; k < 4; k++) {
      for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++) {

          indx=Ind(i+l-1, j+m-1, k+n-1);
          if ((indx >= 0) && (indx <_x*_y*_z)) {
            //P(indx) += b*B(i, s) * B(j, t) * B(k, u)*w;
            P(indx) += b*N(i+l-1,x,_x) * N(j+m-1,y,_y) * N(k+n-1,z,_z)*w;
          } else {
            //cerr<<"out of range:"<< i+l-1<<" "<<j+m-1<<" "<<k+n-1<<endl;
          }
          //cerr<<i+l-1<<" "<<j+m-1<<" "<<k+n-1<<endl;
          for (kk = 0; kk < 4; kk++) {
            for (jj = 0; jj < 4; jj++) {
              for (ii = 0; ii < 4; ii++) {
                indx1=Ind(ii+l-1, jj+m-1, kk+n-1);
                if ((indx >= 0) && (indx <_x*_y*_z) && (indx1 >= 0) && (indx1 <_x*_y*_z)) {
                  M(indx,indx1)
                  //+= B(i, s) * B(j, t) * B(k, u) * B(ii, s) * B(jj, t) * B(kk, u)*w;
                  += N(i+l-1,x,_x) * N(j+m-1,y,_y) * N(k+n-1,z,_z)*N(ii+l-1,x,_x) * N(jj+m-1,y,_y) * N(kk+n-1,z,_z)*w;
                } else {
                  //cerr<<"out of range:"<< i+l-1<<" "<<j+m-1<<" "<<k+n-1<< ii+l-1<<" "<<jj+m-1<<" "<<kk+n-1<<endl;
                }
              }
            }
          }

        }
      }
    }

  }

  for (i=0;i<_x*_y*_z; i++)
    for (j=0;j<_x*_y*_z; j++) {
      T(i)+=M(i,j);
    }
  //M.Print();

  cerr<<"Right side:"<<endl;
  //P.Print();

  cerr<<"Singular value decomposition:"<<endl;

  double **aa, *ww, **vv, *pp, *xx;
  //int *ind;
  int rows = _x*_y*_z;

  // Allocate memory
  aa = dmatrix(1, rows, 1, rows);
  vv = dmatrix(1, rows, 1, rows);
  ww = dvector(1, rows);
  pp = dvector(1, rows);
  xx = dvector(1, rows);

  // Convert matrix to NR format
  M.Matrix2NR(aa);
  P.Vector2NR(pp);

  //ind = ivector(1, rows);

  svdcmp(aa, rows, rows, ww, vv);

 // for (j = 1; j <= rows; j++)  cerr <<ww[j]<<endl;


  svbksb(aa,ww,vv,rows,rows,pp,xx);
  P.NR2Vector(xx);
  // Convert NR format back
  //M.NR2Matrix();

  // Deallocate memory
  free_dmatrix(aa, 1, rows, 1, rows);
  free_dmatrix(vv, 1, rows, 1, rows);
  free_dvector(ww, 1,rows);
  free_dvector(pp, 1,rows);
  free_dvector(xx, 1,rows);
  //free_ivector( ind, 1, rows);
  cerr<<"Resulting  control points:"<<endl;
  for (k = 0; k < _z; k++) {
    for (j = 0; j < _y; j++) {
      for (i = 0; i < _x; i++) {
        cerr<<P(Ind(i,j,k))<<" ";
      }
      cerr<<" ";
    }
    cerr<<endl;
  }

  for (k = 0; k < _z; k++) {
    for (j = 0; j < _y; j++) {
      for (i = 0; i < _x; i++) {
        _data[k][j][i] = P(Ind(i,j,k));
      }
    }
  }
}



void irtkBSplineBiasField::Interpolate(double* dbias)
{
  irtkRealImage coeffs;

  ComputeCoefficients(dbias, coeffs);

  for (int z = 0; z < _z; z++)
    for (int y = 0; y < _y; y++)
      for (int x = 0; x < _x; x++) {
        _data[z][y][x] = coeffs(x, y, z);
      }
}

void irtkBSplineBiasField::Subdivide()
{
  if (_z == 1) {
    this->Subdivide2D();
  } else {
    this->Subdivide3D();
  }
}

void irtkBSplineBiasField::Subdivide2D()
{
  int i, j, i1, j1, i2, j2;

  // Weights for subdivision
  double w[2][3];
  w[1][0] = 0;
  w[1][1] = 1.0/2.0;
  w[1][2] = 1.0/2.0;
  w[0][0] = 1.0/8.0;
  w[0][1] = 6.0/8.0;
  w[0][2] = 1.0/8.0;

  // Allocate memory for new control points
  double ***x = NULL;
  x = this->Allocate(x, 2*_x-1, 2*_y-1, _z);

  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (i1 = 0; i1 < 2; i1++) {
        for (j1 = 0; j1 < 2; j1++) {
          x[0][2*j+j1][2*i+i1] = 0;
          for (i2 = 0; i2 < 3; i2++) {
            for (j2 = 0; j2 < 3; j2++) {
              x[0][2*j+j1][2*i+i1] += w[i1][i2] *
                                      w[j1][j2] * _data[0][j+j2-1][i+i2-1];
            }
          }
        }
      }
    }
  }

  // Deallocate points
  this->Deallocate(_data, _x, _y, _z);

  // Update pointers to control points
  _data = x;

  // Increase number of control points
  _x = 2*_x - 1;
  _y = 2*_y - 1;
  _z = 1;

  // Recalculate control point spacing
  _dx /= 2.0;
  _dy /= 2.0;

  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBSplineBiasField::Subdivide3D()
{
  int i, j, k, i1, j1, k1, i2, j2, k2;

  // Weights for subdivision
  double w[2][3];
  w[1][0] = 0;
  w[1][1] = 1.0/2.0;
  w[1][2] = 1.0/2.0;
  w[0][0] = 1.0/8.0;
  w[0][1] = 6.0/8.0;
  w[0][2] = 1.0/8.0;

  // Allocate memory for new control points
  double ***x = NULL;
  x = this->Allocate(x, 2*_x-1, 2*_y-1, 2*_z-1);

  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        for (i1 = 0; i1 < 2; i1++) {
          for (j1 = 0; j1 < 2; j1++) {
            for (k1 = 0; k1 < 2; k1++) {
              x[2*k+k1][2*j+j1][2*i+i1] = 0;
              for (i2 = 0; i2 < 3; i2++) {
                for (j2 = 0; j2 < 3; j2++) {
                  for (k2 = 0; k2 < 3; k2++) {
                    x[2*k+k1][2*j+j1][2*i+i1] += w[i1][i2] *
                                                 w[j1][j2] * w[k1][k2] * _data[k+k2-1][j+j2-1][i+i2-1];
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
  this->Deallocate(_data, _x, _y, _z);

  // Update pointers to control points
  _data = x;

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
}

double irtkBSplineBiasField::InitialAntiCausalCoefficient(double c[], int DataLength, double z)
{
  /* this initialization corresponds to mirror boundaries */
  return((z / (z * z - 1.0)) * (z * c[DataLength - 2] + c[DataLength - 1]));
}

double irtkBSplineBiasField::InitialCausalCoefficient(double c[], int DataLength, double z, double Tolerance)
{
  double Sum, zn, z2n, iz;
  int n, Horizon;

  /* this initialization corresponds to mirror boundaries */
  Horizon = DataLength;
  if (Tolerance > 0.0) {
    Horizon = (int)ceil(log(Tolerance) / log(fabs(z)));
  }
  if (Horizon < DataLength) {
    /* accelerated loop */
    zn = z;
    Sum = c[0];
    for (n = 1; n < Horizon; n++) {
      Sum += zn * c[n];
      zn *= z;
    }
    return(Sum);
  } else {
    /* full loop */
    zn = z;
    iz = 1.0 / z;
    z2n = pow(z, (double)(DataLength - 1));
    Sum = c[0] + z2n * c[DataLength - 1];
    z2n *= z2n * iz;
    for (n = 1; n <= DataLength - 2; n++) {
      Sum += (zn + z2n) * c[n];
      zn *= z;
      z2n *= iz;
    }
    return(Sum / (1.0 - zn * zn));
  }
}

void irtkBSplineBiasField::ConvertToInterpolationCoefficients(double* c, int DataLength, double* z, int NbPoles, double Tolerance)
{
  double Lambda = 1.0;
  int n, k;

  /* special case required by mirror boundaries */
  if (DataLength == 1) {
    return;
  }

  /* compute the overall gain */
  for (k = 0; k < NbPoles; k++) {
    Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
  }

  /* apply the gain */
  for (n = 0; n < DataLength; n++) {
    c[n] *= Lambda;
  }

  /* loop over all poles */
  for (k = 0; k < NbPoles; k++) {
    /* causal initialization */
    c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
    /* causal recursion */
    for (n = 1; n < DataLength; n++) {
      c[n] += z[k] * c[n - 1];
    }
    /* anticausal initialization */
    c[DataLength - 1] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
    /* anticausal recursion */
    for (n = DataLength - 2; 0 <= n; n--) {
      c[n] = z[k] * (c[n + 1] - c[n]);
    }
  }
}

void irtkBSplineBiasField::ComputeCoefficients(double* dbias, irtkRealImage& coeffs)
{
  int x, y, z, NbPoles = 1;
  double Pole[2];

  Pole[0] = sqrt(3.0) - 2.0;

  // Initialize coefficient images.
  coeffs=irtkRealImage(_x, _y, _z);

  // Convert the displacements into interpolation coefficients for each
  // direction.

  // In-place separable process, along x.
  double* data = new double[_x];
  for (z = 0; z < _z; z++) {
    for (y = 0; y < _y; y++) {
      for (x = 0; x < _x; x++) {
        int index = x + y*_x + z*_x*_y;

        data[x] = dbias[index];
      }

      ConvertToInterpolationCoefficients(data, _x, Pole, NbPoles,
                                         DBL_EPSILON);

      for (x = 0; x < _x; x++) {
        coeffs(x, y, z) = data[x];
      }
    }
  }
  delete[] data;

  // In-place separable process, along y.
  data = new double[_y];
  for (z = 0; z < _z; z++) {
    for (x = 0; x < _x; x++) {
      for (y = 0; y < _y; y++) {
        data[y] = coeffs(x, y, z);
      }

      ConvertToInterpolationCoefficients(data, _y, Pole, NbPoles,
                                         DBL_EPSILON);

      for (y = 0; y < _y; y++) {
        coeffs(x, y, z) = data[y];
      }
    }
  }
  delete[] data;

  // In-place separable process, along z.
  data = new double[_z];
  for (y = 0; y < _y; y++) {
    for (x = 0; x < _x; x++) {
      for (z = 0; z < _z; z++) {
        data[z] = coeffs(x, y, z);
      }

      ConvertToInterpolationCoefficients(data, _z, Pole, NbPoles,
                                         DBL_EPSILON);

      for (z = 0; z < _z; z++) {
        coeffs(x, y, z) = data[z];
      }
    }
  }
  delete[] data;
}

void irtkBSplineBiasField::Read(char *name)
{
  int i, j, k;
  unsigned int magic_no;
  unsigned int trans_type;

  // Open file
  irtkCifstream from;
  from.Open(name);

  // Read magic no. for transformations
  from.ReadAsUInt(&magic_no, 1, 0);

  // Read transformation type
  from.ReadAsUInt(&trans_type, 1);

  if ((magic_no != IRTKBIASFIELD_MAGIC) && (trans_type != IRTKBIASFIELD_BSPLINE)) {
    cerr << "irtkBSplineBiasField::Read: File format not recognized" << endl;
    exit(1);
  }

  // Free memory if necessary
  _data = this->Deallocate(_data, _x, _y, _z);

  // Read no of control points
  from.ReadAsInt(&_x, 1);
  from.ReadAsInt(&_y, 1);
  from.ReadAsInt(&_z, 1);

  // Read orientation
  from.ReadAsDouble(_xaxis, 3);
  from.ReadAsDouble(_yaxis, 3);
  from.ReadAsDouble(_zaxis, 3);

  // Read spacing
  from.ReadAsDouble(&_dx, 1);
  from.ReadAsDouble(&_dy, 1);
  from.ReadAsDouble(&_dz, 1);

  // Read origin
  from.ReadAsDouble(&_origin._x, 1);
  from.ReadAsDouble(&_origin._y, 1);
  from.ReadAsDouble(&_origin._z, 1);

  // Initialize control points
  _data = this->Allocate(_data, _x, _y, _z);

  // Initialize control points
  for (i = -2; i < _x+2; i++) {
    for (j = -2; j < _y+2; j++) {
      for (k = -2; k < _z+2; k++) {
        _data[k][j][i] = 0;
      }
    }
  }

  // Allocate temporary memory
  double *data = new double[_x*_y*_z];

  // Read control point data
  from.ReadAsDouble(data, _x*_y*_z, 512);

  // Convert control point data
  int index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        _data[k][j][i] = data[index];
        index++;
      }
    }
  }

  // Delete memory
  delete []data;

  // Close file stream
  from.Close();

  // Update transformation matrix
  this->UpdateMatrix();
}


void irtkBSplineBiasField::Write(char *name)
{
  int i, j, k;

  // Open file
  irtkCofstream to;
  to.Open(name);

  // Write magic no. for transformations
  unsigned int magic_no = IRTKBIASFIELD_MAGIC;
  to.WriteAsUInt(&magic_no, 1, 0);

  // Write transformation type
  unsigned int trans_type = IRTKBIASFIELD_BSPLINE;
  to.WriteAsUInt(&trans_type, 1);

  // Write no of control points
  to.WriteAsInt(&_x, 1);
  to.WriteAsInt(&_y, 1);
  to.WriteAsInt(&_z, 1);

  // Write orientation
  to.WriteAsDouble(_xaxis, 3);
  to.WriteAsDouble(_yaxis, 3);
  to.WriteAsDouble(_zaxis, 3);

  // Write spacing
  to.WriteAsDouble(&_dx, 1);
  to.WriteAsDouble(&_dy, 1);
  to.WriteAsDouble(&_dz, 1);

  // Write origin
  to.WriteAsDouble(&_origin._x, 1);
  to.WriteAsDouble(&_origin._y, 1);
  to.WriteAsDouble(&_origin._z, 1);

  // Write zeros
  for (i = 140; i < 512; i++) {
    // Write 0
    unsigned char zero = 0;
    to.WriteAsUChar(&zero, 1, i);
  }

  // Allocate temporary memory
  double *data = new double[_x*_y*_z];

  // Convert control point data
  int index = 0;
  for (i = 0; i < _x; i++) {
    for (j = 0; j < _y; j++) {
      for (k = 0; k < _z; k++) {
        data[index]  = _data[k][j][i];
        index++;
      }
    }
  }

  // Write control point data
  to.WriteAsDouble(data, _x*_y*_z, 512);

  // Delete memory
  delete []data;

  // Close file stream
  to.Close();
}

void irtkBSplineBiasField::Print()
{
  cerr<<endl<<"BSplineBiasField info:" << endl;
  // Write no. of control points
  cout << "Control points: " << _x << " x " << _y << " x " << _z << endl;
  cout << "Spacing: " << _dx << " x " << _dy << " x " << _dz << endl;
  cout << "Origin: " << _origin._x << " " << _origin._y << " " << _origin._z << " " << endl;
  cout << "Orientation: " << _xaxis[0] << " " << _xaxis[1] << " "
       << _xaxis[2] << " " << _yaxis[0] << " " << _yaxis[1] << " "
       << _yaxis[2] << " " << _zaxis[0] << " " << _zaxis[1] << " " << _zaxis[2] << endl;
  cout << "B spline control points:" << endl;

  int i,j,k;
  for (i=0; i<_z; i++) {
    for (j=0; j<_y; j++) {
      for (k=01; k<_x; k++) cerr << _data[i][j][k]<<" ";
      cerr<<"     ";
    }
    cerr<<endl;
  }
  cerr<<endl;
}

double irtkBSplineBiasField::N(int i, double u, int L)
{
	if((i<0)||(i>L-1)) return 0;
	if((u<0)||(u>L-1)) return 0;
	int l= (int) floor(u);
	if(l==L-1) l=L-2;
	double t = u-l;
    double t3=t*t*t, t2=t*t;

    if((l==0)&&(L>3))
	{
    	switch (i-l+1) {
    	case 0:
    		return 0;
		case 1:
			return (9*t3 - 18*t2 + 12)/12.0;
		case 2:
		    return (-11*t3 + 18*t2 )/12.0;
		case 3:
		    return (t3)/6.0;
		}
		return 0;
	}

	  if((l==L-2)&&(L>3))
	  {
	    switch (i-l+1) {
	    case 0:
	      return (-2*t3 + 6*t2 -6*t + 2)/12.0;;
	    case 1:
	      return (11*t3 - 15*t2 -3*t + 7)/12.0;
	    case 2:
	      return (-9*t3 + 9*t2 + 9*t + 3)/12.0;
	    case 3:
	      return 0;
	    }
	    return 0;
	  }

	 if((l==1)&&(L>4))
	 {
	    switch (i-l+1) {
	    case 0:
	      return (-3*t3 + 9*t2 -9*t +3)/12.0;
	    case 1:
	      return (7*t3 - 15*t2 + 3*t +7)/12.0;
	    case 2:
	      return (-6*t3 + 6*t2 + 6*t + 2 )/12.0;
	    case 3:
	      return (t3)/6.0;
	    }
	    return 0;
	  }

	  if((l==L-3)&&(L>4))
	  {
	    switch (i-l+1) {
	    case 0:
	      return (-2*t3 + 6*t2 -6*t + 2)/12.0;;
	    case 1:
	      return (6*t3 - 12*t2  + 8)/12.0;
	    case 2:
	      return (-7*t3 + 6*t2 + 6*t + 2)/12.0;
	    case 3:
	      return (3*t3)/12.0;
	    }
	    return 0;
	  }
	  if((l==1)&&(L==4))
	  {
	    switch (i-l+1) {
	    case 0:
	      return (-3*t3 + 9*t2 -9*t + 3)/12.0;
	    case 1:
	      return (7*t3 - 15*t2  +3*t + 7)/12.0;
	    case 2:
	      return (-7*t3 + 6*t2 + 6*t + 2)/12.0;
	    case 3:
	      return (3*t3)/12.0;
	    }
	    return 0;
	  }

	  if((l==0)&&(L==3))
	  {
	    switch (i-l+1) {
	    case 0:
	      return 0;
	    case 1:
	      return (3*t3 - 6*t2  + 4)/4.0;
	    case 2:
	      return (-4*t3 + 6*t2 )/4.0;
	    case 3:
	      return (t3)/4.0;
	    }
	    return 0;
	  }

	  if((l==1)&&(L==3))
	  {
	    switch (i-l+1) {
	    case 0:
	      return (-1*t3 + 3*t2 - 3*t + 1)/4.0;
	    case 1:
	      return (4*t3 - 6*t2  + 2)/4.0;
	    case 2:
	      return (-3*t3 + 3*t2 + 3*t + 1)/4.0;
	    case 3:
	      return 0;
	    }
	    return 0;
	  }

	  if((l==0)&&(L==2))
	  {
	    switch (i-l+1) {
	    case 0:
	      return 0;
	    case 1:
	      return 2*t3 - 3*t2  + 1;
	    case 2:
	      return -2*t3 + 3*t2;
	    case 3:
	      return 0;
	    }
	    return 0;
	  }

	  if((l<L-3)&&(l>1))
	  {
	    switch (i-l+1) {
	    case 0:
	      return (-t3+3*t2-3*t+1)/6.0;
	    case 1:
	      return (3*t3 - 6*t2 + 4)/6.0;
	    case 2:
	      return (-3*t3 + 3*t2 + 3*t + 1)/6.0;
	    case 3:
	      return (t3)/6.0;
	    }
	    return 0;
	  }
	  return 0;
  }

