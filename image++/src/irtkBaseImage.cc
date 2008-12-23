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

#include <irtkNIFTI.h>

irtkBaseImage::irtkBaseImage()
{
  _x  = 0;
  _y  = 0;
  _z  = 0;
  _t  = 0;

  // Default voxel size
  _dx = 1;
  _dy = 1;
  _dz = 1;
  _dt = 1;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  _matI2W = irtkMatrix(4, 4);
  _matW2I = irtkMatrix(4, 4);

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Update transformation matrix
  this->UpdateMatrix();
}

irtkBaseImage::irtkBaseImage(int x, int y, int z, int t)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = t;

  // Default voxel size
  _dx = 1;
  _dy = 1;
  _dz = 1;
  _dt = 1;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  _matI2W = irtkMatrix(4, 4);
  _matW2I = irtkMatrix(4, 4);

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Update transformation matrix
  this->UpdateMatrix();
}

irtkBaseImage::irtkBaseImage(int x, int y, int z, double dx, double dy, double dz)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = 1;

  // Default voxel size
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = 1;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  _matI2W = irtkMatrix(4, 4);
  _matW2I = irtkMatrix(4, 4);

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Update transformation matrix
  this->UpdateMatrix();
}

irtkBaseImage::irtkBaseImage(int x, int y, int z, double dx, double dy, double dz,
                             double *xaxis, double *yaxis, double *zaxis)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _z  = 1;

  // Default voxel size
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dz = 1;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  _matI2W = irtkMatrix(4, 4);
  _matW2I = irtkMatrix(4, 4);

  // Default x-axis
  _xaxis[0] = xaxis[0];
  _xaxis[1] = xaxis[1];
  _xaxis[2] = xaxis[2];

  // Default y-axis
  _yaxis[0] = yaxis[0];
  _yaxis[1] = yaxis[1];
  _yaxis[2] = yaxis[2];

  // Default z-axis
  if (zaxis != NULL) {
    _zaxis[0] = zaxis[0];
    _zaxis[1] = zaxis[1];
    _zaxis[2] = zaxis[2];
  } else {
    // Update z-axis
    _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
    _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
    _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

irtkBaseImage::irtkBaseImage(int x, int y, int z, double dx, double dy, double dz, irtkPoint origin,
                             double *xaxis, double *yaxis, double *zaxis)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = 1;

  // Default voxel size
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = 1;

  // Default origin
  _origin = origin;
  _torigin = 0;

  _matI2W = irtkMatrix(4, 4);
  _matW2I = irtkMatrix(4, 4);

  // Default x-axis
  _xaxis[0] = xaxis[0];
  _xaxis[1] = xaxis[1];
  _xaxis[2] = xaxis[2];

  // Default y-axis
  _yaxis[0] = yaxis[0];
  _yaxis[1] = yaxis[1];
  _yaxis[2] = yaxis[2];

  // Default z-axis
  if (zaxis != NULL) {
    _zaxis[0] = zaxis[0];
    _zaxis[1] = zaxis[1];
    _zaxis[2] = zaxis[2];
  } else {
    // Update z-axis
    _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
    _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
    _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

irtkBaseImage::irtkBaseImage(const irtkBaseImage &image)
{
  // Copy image dimensions
  _x = image._x;
  _y = image._y;
  _z = image._z;
  _t = image._t;

  // Copy voxel dimensions
  _dx = image._dx;
  _dy = image._dy;
  _dz = image._dz;
  _dt = image._dt;

  // Copy origin
  _origin = image._origin;
  _torigin = image._torigin;
  _matI2W = image._matI2W;
  _matW2I = image._matW2I;

  // Copy x-axis
  _xaxis[0] = image._xaxis[0];
  _xaxis[1] = image._xaxis[1];
  _xaxis[2] = image._xaxis[2];

  // Copy y-axis
  _yaxis[0] = image._yaxis[0];
  _yaxis[1] = image._yaxis[1];
  _yaxis[2] = image._yaxis[2];

  // Copy z-axis
  _zaxis[0] = image._zaxis[0];
  _zaxis[1] = image._zaxis[1];
  _zaxis[2] = image._zaxis[2];
}

irtkBaseImage::~irtkBaseImage()
{
  _x  = 0;
  _y  = 0;
  _z  = 0;
  _t  = 0;

  // Default voxel size
  _dx = 1;
  _dy = 1;
  _dz = 1;
  _dt = 1;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  _matI2W = irtkMatrix(4, 4);
  _matW2I = irtkMatrix(4, 4);

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBaseImage::UpdateMatrix()
{
  // Note that the update of zaxis is taken place in the
  // Initialize routines, in case a specific zaxis is passed.

  // Update image to world coordinate system matrix
  _matI2W.Ident();
  _matI2W(0, 0) = _xaxis[0];
  _matI2W(1, 0) = _xaxis[1];
  _matI2W(2, 0) = _xaxis[2];
  _matI2W(0, 1) = _yaxis[0];
  _matI2W(1, 1) = _yaxis[1];
  _matI2W(2, 1) = _yaxis[2];
  _matI2W(0, 2) = _zaxis[0];
  _matI2W(1, 2) = _zaxis[1];
  _matI2W(2, 2) = _zaxis[2];

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

  _matI2W = tmp3 * (_matI2W * (tmp2 * tmp1));

  // Update world to image coordinate system matrix
  _matW2I.Ident();
  _matW2I(0, 0) = _xaxis[0];
  _matW2I(0, 1) = _xaxis[1];
  _matW2I(0, 2) = _xaxis[2];
  _matW2I(1, 0) = _yaxis[0];
  _matW2I(1, 1) = _yaxis[1];
  _matW2I(1, 2) = _yaxis[2];
  _matW2I(2, 0) = _zaxis[0];
  _matW2I(2, 1) = _zaxis[1];
  _matW2I(2, 2) = _zaxis[2];

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

  _matW2I = tmp1 * (tmp2 * (_matW2I * tmp3));
}

void irtkBaseImage::Initialize(int x, int y, int z, int t)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = t;

  // Default voxel size
  _dx = 1;
  _dy = 1;
  _dz = 1;
  _dt = 1;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBaseImage::Initialize(int x, int y, int z, double dx, double dy, double dz)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = 1;
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = 1;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBaseImage::Initialize(int x, int y, int z, int t, double dx, double dy, double dz, double dt)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = t;
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = dt;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  // Default x-axis
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;

  // Default y-axis
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;

  // Default z-axis
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBaseImage::Initialize(int x, int y, int z, double dx, double dy, double dz, const double *xaxis, const double *yaxis, const double *zaxis)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = 1;
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = 1;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  // Default x-axis
  _xaxis[0] = xaxis[0];
  _xaxis[1] = xaxis[1];
  _xaxis[2] = xaxis[2];

  // Default y-axis
  _yaxis[0] = yaxis[0];
  _yaxis[1] = yaxis[1];
  _yaxis[2] = yaxis[2];

  // Default z-axis
  if (zaxis != NULL) {
    _zaxis[0] = zaxis[0];
    _zaxis[1] = zaxis[1];
    _zaxis[2] = zaxis[2];
  } else {
    // Update z-axis
    _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
    _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
    _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBaseImage::Initialize(int x, int y, int z, int t, double dx, double dy, double dz, double dt, const double *xaxis, const double *yaxis, const double *zaxis)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = t;
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = dt;

  // Default origin
  _origin = irtkPoint();
  _torigin = 0;

  // Default x-axis
  _xaxis[0] = xaxis[0];
  _xaxis[1] = xaxis[1];
  _xaxis[2] = xaxis[2];

  // Default y-axis
  _yaxis[0] = yaxis[0];
  _yaxis[1] = yaxis[1];
  _yaxis[2] = yaxis[2];

  // Default z-axis
  if (zaxis != NULL) {
    _zaxis[0] = zaxis[0];
    _zaxis[1] = zaxis[1];
    _zaxis[2] = zaxis[2];
  } else {
    // Update z-axis
    _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
    _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
    _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBaseImage::Initialize(int x, int y, int z, double dx, double dy, double dz, irtkPoint origin,
                               const double *xaxis, const double *yaxis, const double *zaxis)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = 1;
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = 1;

  // Default origin
  _origin = origin;
  _torigin = 0;

  // Default x-axis
  _xaxis[0] = xaxis[0];
  _xaxis[1] = xaxis[1];
  _xaxis[2] = xaxis[2];

  // Default y-axis
  _yaxis[0] = yaxis[0];
  _yaxis[1] = yaxis[1];
  _yaxis[2] = yaxis[2];

  // Default z-axis
  if (zaxis != NULL) {
    _zaxis[0] = zaxis[0];
    _zaxis[1] = zaxis[1];
    _zaxis[2] = zaxis[2];
  } else {
    // Update z-axis
    _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
    _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
    _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBaseImage::Initialize(int x, int y, int z, int t, double dx, double dy, double dz, double dt, irtkPoint origin, double torigin, const double *xaxis, const double *yaxis, const double *zaxis)
{
  _x  = x;
  _y  = y;
  _z  = z;
  _t  = t;
  _dx = dx;
  _dy = dy;
  _dz = dz;
  _dt = dt;

  // Default origin
  _origin = origin;
  _torigin = torigin;

  // Default x-axis
  _xaxis[0] = xaxis[0];
  _xaxis[1] = xaxis[1];
  _xaxis[2] = xaxis[2];

  // Default y-axis
  _yaxis[0] = yaxis[0];
  _yaxis[1] = yaxis[1];
  _yaxis[2] = yaxis[2];

  // Default z-axis
  if (zaxis != NULL) {
    _zaxis[0] = zaxis[0];
    _zaxis[1] = zaxis[1];
    _zaxis[2] = zaxis[2];
  } else {
    // Update z-axis
    _zaxis[0] = _xaxis[1]*_yaxis[2] - _xaxis[2]*_yaxis[1];
    _zaxis[1] = _xaxis[2]*_yaxis[0] - _xaxis[0]*_yaxis[2];
    _zaxis[2] = _xaxis[0]*_yaxis[1] - _xaxis[1]*_yaxis[0];
  }

  // Update transformation matrix
  this->UpdateMatrix();
}

void irtkBaseImage::Initialize(const irtkBaseImage &image)
{
  // Copy image dimensions
  _x = image._x;
  _y = image._y;
  _z = image._z;
  _t = image._t;

  // Copy voxel dimensions
  _dx = image._dx;
  _dy = image._dy;
  _dz = image._dz;
  _dt = image._dt;

  // Copy origin and matrix
  _origin = image._origin;
  _matI2W = image._matI2W;
  _matW2I = image._matW2I;

  // Copy axes
  _xaxis[0] = image._xaxis[0];
  _xaxis[1] = image._xaxis[1];
  _xaxis[2] = image._xaxis[2];
  _yaxis[0] = image._yaxis[0];
  _yaxis[1] = image._yaxis[1];
  _yaxis[2] = image._yaxis[2];
  _zaxis[0] = image._zaxis[0];
  _zaxis[1] = image._zaxis[1];
  _zaxis[2] = image._zaxis[2];
}

irtkBaseImage &irtkBaseImage::operator=(const irtkBaseImage &image)
{
  this->Initialize(image);
  return *this;
}

Bool irtkBaseImage::operator==(const irtkBaseImage &image)
{
  return ((_x  == image._x)  && (_y  == image._y)  && (_z  == image._z) && (_t  == image._t) &&
          (_dx == image._dx) && (_dy == image._dy) && (_dz == image._dz) && (_dt == image._dt) &&
          (_origin == image._origin) && (_torigin == image._torigin));
}

void irtkBaseImage::Orientation(int &i, int &j, int &k) const
{

#ifdef HAS_NIFTI

  // Work out orientation of axis
  mat44 mat_44;
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      mat_44.m[i][j] = this->GetImageToWorldMatrix()(i, j);
    }
  }
  nifti_mat44_to_orientation(mat_44, &i, &j, &k);

#else

  cerr << "irtkBaseImage::Orientation: Requires NIFTI support. Please recompile with NIFTI enabled" << endl;
  exit(1);

#endif
}

void irtkBaseImage::Print()
{
  // Print image dimensions
  cout << "Image size is " << _x << " " << _y << " " << _z << " " << _t << endl;
  // Print voxel dimensions
  cout << "Voxel size is " << _dx << " " << _dy << " "
       << _dz << " " << _dt << endl;
  // Print origin
  cout << "Image origin is " << _origin << " " << _torigin << endl;
  // Print x-axis
  cout << "X-axis is " << _xaxis[0] << " " << _xaxis[1] << " " << _xaxis[2] << endl;
  // Print x-axis
  cout << "Y-axis is " << _yaxis[0] << " " << _yaxis[1] << " " << _yaxis[2] << endl;
  // Print x-axis
  cout << "Z-axis is " << _zaxis[0] << " " << _zaxis[1] << " " << _zaxis[2] << endl;
}
