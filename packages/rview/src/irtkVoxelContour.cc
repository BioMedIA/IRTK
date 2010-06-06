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
#include <irtkTransformation.h>
#include <irtkRegistration.h>

#ifdef __APPLE__
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#ifndef WIN32
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <stack>

#include <irtkRView.h>

typedef struct
{
  int x;
  int y;
  int z;
}
irtkLocation;

irtkVoxelContour::irtkVoxelContour()
{
  _raster = new irtkGreyImage;
  _rview  = NULL;
  _currentSize = 0;
  _totalSize = 0;
  _current = 0;
}

void irtkVoxelContour::Initialise(irtkRView *rview, irtkGreyImage* viewer)
{
  int i, i1, j1, k1, i2, j2, k2, axis[3];
  double x1, y1, z1, x2, y2, z2;

  // Pointer to rview
  _rview = rview;

  // Figure out which plane we are in
  viewer->Orientation(i1, j1, k1);
  _rview->_targetImage->Orientation(i2, j2, k2);
  axis[0] = 0;
  axis[1] = 1;
  axis[2] = 2;
  switch (i1) {
    case IRTK_L2R:
    case IRTK_R2L:
      if ((i2 == IRTK_L2R) || (i2 == IRTK_R2L)) {
        axis[0] = 0;
      }
      if ((j2 == IRTK_L2R) || (j2 == IRTK_R2L)) {
        axis[0] = 1;
      }
      if ((k2 == IRTK_L2R) || (k2 == IRTK_R2L)) {
        axis[0] = 2;
      }
      break;
    case IRTK_A2P:
    case IRTK_P2A:
      if ((i2 == IRTK_A2P) || (i2 == IRTK_P2A)) {
        axis[0] = 0;
      }
      if ((j2 == IRTK_A2P) || (j2 == IRTK_P2A)) {
        axis[0] = 1;
      }
      if ((k2 == IRTK_A2P) || (k2 == IRTK_P2A)) {
        axis[0] = 2;
      }
      break;
    case IRTK_S2I:
    case IRTK_I2S:
      if ((i2 == IRTK_S2I) || (i2 == IRTK_I2S)) {
        axis[0] = 0;
      }
      if ((j2 == IRTK_S2I) || (j2 == IRTK_I2S)) {
        axis[0] = 1;
      }
      if ((k2 == IRTK_S2I) || (k2 == IRTK_I2S)) {
        axis[0] = 2;
      }
      break;
  }
  switch (j1) {
    case IRTK_L2R:
    case IRTK_R2L:
      if ((i2 == IRTK_L2R) || (i2 == IRTK_R2L)) {
        axis[1] = 0;
      }
      if ((j2 == IRTK_L2R) || (j2 == IRTK_R2L)) {
        axis[1] = 1;
      }
      if ((k2 == IRTK_L2R) || (k2 == IRTK_R2L)) {
        axis[1] = 2;
      }
      break;
    case IRTK_A2P:
    case IRTK_P2A:
      if ((i2 == IRTK_A2P) || (i2 == IRTK_P2A)) {
        axis[1] = 0;
      }
      if ((j2 == IRTK_A2P) || (j2 == IRTK_P2A)) {
        axis[1] = 1;
      }
      if ((k2 == IRTK_A2P) || (k2 == IRTK_P2A)) {
        axis[1] = 2;
      }
      break;
    case IRTK_S2I:
    case IRTK_I2S:
      if ((i2 == IRTK_S2I) || (i2 == IRTK_I2S)) {
        axis[1] = 0;
      }
      if ((j2 == IRTK_S2I) || (j2 == IRTK_I2S)) {
        axis[1] = 1;
      }
      if ((k2 == IRTK_S2I) || (k2 == IRTK_I2S)) {
        axis[1] = 2;
      }
      break;
  }
  switch (k1) {
    case IRTK_L2R:
    case IRTK_R2L:
      if ((i2 == IRTK_L2R) || (i2 == IRTK_R2L)) {
        axis[2] = 0;
      }
      if ((j2 == IRTK_L2R) || (j2 == IRTK_R2L)) {
        axis[2] = 1;
      }
      if ((k2 == IRTK_L2R) || (k2 == IRTK_R2L)) {
        axis[2] = 2;
      }
      break;
    case IRTK_A2P:
    case IRTK_P2A:
      if ((i2 == IRTK_A2P) || (i2 == IRTK_P2A)) {
        axis[2] = 0;
      }
      if ((j2 == IRTK_A2P) || (j2 == IRTK_P2A)) {
        axis[2] = 1;
      }
      if ((k2 == IRTK_A2P) || (k2 == IRTK_P2A)) {
        axis[2] = 2;
      }
      break;
    case IRTK_S2I:
    case IRTK_I2S:
      if ((i2 == IRTK_S2I) || (i2 == IRTK_I2S)) {
        axis[2] = 0;
      }
      if ((j2 == IRTK_S2I) || (j2 == IRTK_I2S)) {
        axis[2] = 1;
      }
      if ((k2 == IRTK_S2I) || (k2 == IRTK_I2S)) {
        axis[2] = 2;
      }
      break;
  }
  // Get attributes of viewer
  irtkImageAttributes attr = _rview->_targetImage->GetImageAttributes();

  // Compute position of origin in target voxel coordinates
  x1 = attr._xorigin;
  y1 = attr._yorigin;
  z1 = attr._zorigin;
  _rview->_targetImage->WorldToImage(x1, y1, z1);
  x2 = attr._xorigin;
  y2 = attr._yorigin;
  z2 = attr._zorigin;
  _rview->_targetImage->WorldToImage(x2, y2, z2);

  // Create image
  if ((axis[0] == 0) && (axis[1] == 1)) {
    _raster->Initialize(attr);
    x1 = x2;
    y1 = y2;
    z1 = z2;
  }
  if ((axis[0] == 1) && (axis[1] == 0)) {
    irtkImageAttributes tmp;
    tmp._y = attr._x;
    tmp._x = attr._y;
    tmp._z = attr._z;
    tmp._dy = attr._dx;
    tmp._dx = attr._dy;
    tmp._dz = attr._dz;
    for (i = 0; i < 3; i++) {
      tmp._yaxis[i] = attr._xaxis[i];
      tmp._xaxis[i] = attr._yaxis[i];
      tmp._zaxis[i] = attr._zaxis[i];
    }
    _raster->Initialize(tmp);
    x1 = y2;
    y1 = x2;
    z1 = z2;
  }
  if ((axis[0] == 0) && (axis[1] == 2)) {
    irtkImageAttributes tmp;
    tmp._x = attr._x;
    tmp._z = attr._y;
    tmp._y = attr._z;
    tmp._dx = attr._dx;
    tmp._dz = attr._dy;
    tmp._dy = attr._dz;
    for (i = 0; i < 3; i++) {
      tmp._xaxis[i] = attr._xaxis[i];
      tmp._zaxis[i] = attr._yaxis[i];
      tmp._yaxis[i] = attr._zaxis[i];
    }
    _raster->Initialize(tmp);
    x1 = x2;
    y1 = z2;
    z1 = y2;
  }
  if ((axis[0] == 2) && (axis[1] == 0)) {
    irtkImageAttributes tmp;
    tmp._z = attr._x;
    tmp._x = attr._y;
    tmp._y = attr._z;
    tmp._dz = attr._dx;
    tmp._dx = attr._dy;
    tmp._dy = attr._dz;
    for (i = 0; i < 3; i++) {
      tmp._zaxis[i] = attr._xaxis[i];
      tmp._xaxis[i] = attr._yaxis[i];
      tmp._yaxis[i] = attr._zaxis[i];
    }
    _raster->Initialize(tmp);
    x1 = z2;
    y1 = x2;
    z1 = y2;
  }
  if ((axis[0] == 1) && (axis[1] == 2)) {
    irtkImageAttributes tmp;
    tmp._x = attr._y;
    tmp._y = attr._z;
    tmp._z = attr._x;
    tmp._dx = attr._dy;
    tmp._dy = attr._dz;
    tmp._dz = attr._dx;
    for (i = 0; i < 3; i++) {
      tmp._xaxis[i] = attr._yaxis[i];
      tmp._yaxis[i] = attr._zaxis[i];
      tmp._zaxis[i] = attr._xaxis[i];
     }
    _raster->Initialize(tmp);
    x1 = y2;
    y1 = z2;
    z1 = x2;
  }
  if ((axis[0] == 2) && (axis[1] == 1)) {
    irtkImageAttributes tmp;
    tmp._z = attr._x;
    tmp._y = attr._y;
    tmp._x = attr._z;
    tmp._dz = attr._dx;
    tmp._dy = attr._dy;
    tmp._dx = attr._dz;
    for (i = 0; i < 3; i++) {
      tmp._zaxis[i] = attr._xaxis[i];
      tmp._yaxis[i] = attr._yaxis[i];
      tmp._xaxis[i] = attr._zaxis[i];
    }
    _raster->Initialize(tmp);
    x1 = z2;
    y1 = y2;
    z1 = x2;
  }
  _raster->ImageToWorld(x1, y1, z1);
  _rview->_targetImage->ImageToWorld(x2, y2, z2);
  _raster->PutOrigin(x2 - x1, y2 - y1, z2 - z1);

}

void irtkVoxelContour::AddPoint(irtkPoint p, int width)
{
  int x, y, z;

  // Store current paintbrush width
  _width = width;

  // Check if contour has been initialised
  if (_raster == NULL) {
    cerr << "Please, always initialise irtkVoxelContour before adding any points!" << endl;
    exit(1);
  }
  _raster->WorldToImage(p);
  x = round(p._x);
  y = round(p._y);
  z = round(p._z);

  // Check if inside image
  if ((x < 0) || (y < 0) || (z < 0) || (x >= _raster->GetX()) || (y >= _raster->GetY()) || (z >= _raster->GetZ())) return;

  // Check if this is the first point
  if (_currentSize > 0) {
    lineBresenham(_lastx, _lasty, _lastz, x, y, z);
  } else {
    AddPoint(x, y, z);
    _firstx = x;
    _firsty = y;
    _firstz = z;
  }

  _lastx = x;
  _lasty = y;
  _lastz = z;
}


void irtkVoxelContour::AddPoint(int x, int y, int z)
{
  int i, j, k;

  k = (_width-1)/2;
  for (i= -k; i <= k; i++) {
    for (j = -k; j <=k; j++) {
      if ((x+i >= 0) && (y+j >= 0) && (x+i < _raster->GetX()) && (y+j < _raster->GetY())) {
        if (_raster->Get(x+i, y+j, z) == 0) {
          _raster->Put(x+i, y+j, z, _current);
          _currentSize++;
        }
      }
    }
  }
}

void irtkVoxelContour::AddPointSet()
{
  _totalSize += _currentSize;
  _current++;
  _currentSize = 0;
}

void irtkVoxelContour::AddPointSet(irtkPoint p, int width)
{
  this->AddPointSet();
  this->AddPoint(p, width);
}


void irtkVoxelContour::Close(irtkPoint p, int width)
{
  this->AddPoint(p, width);
  lineBresenham(_firstx, _firsty, _firstz, _lastx, _lasty, _lastz);
}

void irtkVoxelContour::FillArea(irtkPoint p)
{
  int x, y, z;

  // Add the filled area as a new point set
  this->AddPointSet();

  // Make sure the point set has at least one point
  if (this->Size() == 0) this->AddPoint(p, 1);

  _raster->WorldToImage(p);
  x = round(p._x);
  y = round(p._y);
  z = round(p._z);

  p._x = x;
  p._y = y;
  p._z = z;
  _raster->ImageToWorld(p._x, p._y, p._z);
  _rview->_targetImage->WorldToImage(p._x, p._y, p._z);

  if ((round(p._x) < 0) || (round(p._x) >= _rview->_targetImage->GetX()) ||
      (round(p._y) < 0) || (round(p._y) >= _rview->_targetImage->GetY()) ||
      (round(p._z) < 0) || (round(p._z) >= _rview->_targetImage->GetZ())) {
    cerr << "Target coordinates out of range: " << x << " " << y << " " << z << " " << endl;
    return;
  }

  // Start recursion
  Fill(x, y, z);
}

void irtkVoxelContour::RegionGrowing(irtkPoint p, int thresholdMin, int thresholdMax, irtkRegionGrowingMode mode)
{
  int x, y, z;

  // Add the filled area as a new point set
  this->AddPointSet();

  // Make sure the point set has at least one point
  if (this->Size() == 0) this->AddPoint(p, 1);

  _raster->WorldToImage(p);
  x = round(p._x);
  y = round(p._y);
  z = round(p._z);

  p._x = x;
  p._y = y;
  p._z = z;
  _raster->ImageToWorld(p._x, p._y, p._z);
  _rview->_targetImage->WorldToImage(p._x, p._y, p._z);

  if ((round(p._x) < 0) || (round(p._x) >= _rview->_targetImage->GetX()) ||
      (round(p._y) < 0) || (round(p._y) >= _rview->_targetImage->GetY()) ||
      (round(p._z) < 0) || (round(p._z) >= _rview->_targetImage->GetZ())) {
    cerr << "Target coordinates out of range: " << x << " " << y << " " << z << " " << endl;
    return;
  }

  // Start region growing
  if (mode == ::RegionGrowing3D) {
    RegionGrowing3D(x, y, z, thresholdMin, thresholdMax);
  } else {
    RegionGrowing2D(x, y, z, thresholdMin, thresholdMax);
  }
}

void irtkVoxelContour::Undo()
{
  int i;
  irtkGreyPixel *ptr;

  if (_current == 0) return;

  ptr = _raster->GetPointerToVoxels();
  for (i = 0; i < _raster->GetNumberOfVoxels(); i++) {
    if (*ptr == _current) *ptr = 0;
    ptr++;
  }
  _currentSize = 0;
  _current--;
}

void irtkVoxelContour::Clear()
{
  int i;
  irtkGreyPixel *ptr;

  ptr = _raster->GetPointerToVoxels();
  for (i = 0; i < _raster->GetNumberOfVoxels(); i++) {
    *ptr = 0;
    ptr++;
  }
  _currentSize = 0;
  _totalSize = 0;
  _current = 0;
}

void irtkVoxelContour::lineBresenham(int x0, int y0, int z0, int x1, int y1, int z1)
{

  int dy = y1 - y0;
  int dx = x1 - x0;
  int stepx, stepy;

  if (z1 != z0) {
    cerr << "lineBresenham: Should not happen" << endl;
  }

  if ((dx == 0) && (dy == 0)) return;

  if (dy < 0) {
    dy = -dy; stepy = -1;
  } else {
    stepy = 1;
  }
  if (dx < 0) {
    dx = -dx; stepx = -1;
  } else {
    stepx = 1;
  }
  dy <<= 1; // dy is now 2*dy
  dx <<= 1; // dx is now 2*dx

  this->AddPoint(x0, y0, z0);

  if (dx > dy) {
    int fraction = dy - (dx >> 1);  // same as 2*dy - dx
    while (x0 != x1) {
      if (fraction >= 0) {
        y0 += stepy;
        fraction -= dx;     // same as fraction -= 2*dx
      }
      x0 += stepx;
      fraction += dy;         // same as fraction -= 2*dy
      AddPoint(x0, y0, z0);
    }
  } else {
    int fraction = dx - (dy >> 1);
    while (y0 != y1) {
      if (fraction >= 0) {
        x0 += stepx;
        fraction -= dy;
      }
      y0 += stepy;
      fraction += dx;
      AddPoint(x0, y0, z0);
    }
  }
}

bool inline irtkVoxelContour::RegionGrowingCriteria(int i, int j, int k, double lowT, double highT)
{
  double x, y, z;

  x = i;
  y = j;
  z = k;
  _raster->ImageToWorld(x, y, z);
  _rview->_targetImage->WorldToImage(x, y, z);
  i = round(x);
  j = round(y);
  k = round(z);
  return ((_rview->_targetImage->GetAsDouble(i, j, k) >= lowT) && (_rview->_targetImage->GetAsDouble(i, j, k) <= highT));
}

void irtkVoxelContour::Fill(int seedX, int seedY, int seedZ)
{
  int x, y, z;
  irtkLocation location;

  // Create stack
  stack<irtkLocation> point_stack;

  // Push seed location on stack
  location.x = seedX;
  location.y = seedY;
  location.z = seedZ;
  point_stack.push(location);
  _raster->Put(location.x, location.y, location.z, _current);

  while (!point_stack.empty()) {
    location = point_stack.top();
    point_stack.pop();
    x = location.x;
    y = location.y;
    z = location.z;
    if (x-1 >= 0) {
      if (_raster->Get(x-1, y, z) == 0) {
        location.x = x-1;
        location.y = y;
        location.z = z;
        point_stack.push(location);
        _raster->Put(x-1, y, z, _current);
        _currentSize++;
      }
    }
    if (x+1 < _raster->GetX()) {
      if (_raster->Get(x+1, y, z) == 0) {
        location.x = x+1;
        location.y = y;
        location.z = z;
        point_stack.push(location);
        _raster->Put(x+1, y, z, _current);
        _currentSize++;
      }
    }
    if (y-1 >= 0) {
      if (_raster->Get(x, y-1, z) == 0) {
        location.x = x;
        location.y = y-1;
        location.z = z;
        point_stack.push(location);
        _raster->Put(x, y-1, z, _current);
        _currentSize++;
      }
    }
    if (y+1 < _raster->GetY()) {
      if (_raster->Get(x, y+1, z) == 0) {
        location.x = x;
        location.y = y+1;
        location.z = z;
        point_stack.push(location);
        _raster->Put(x, y+1, z, _current);
        _currentSize++;
      }
    }
  }
}

void irtkVoxelContour::RegionGrowing2D(int seedX, int seedY, int seedZ, double lowT, double highT)
{
  int x, y, z;
  irtkLocation location;

  // Create stack
  stack<irtkLocation> point_stack;

  // Push seed location on stack
  location.x = seedX;
  location.y = seedY;
  location.z = seedZ;
  point_stack.push(location);
  _raster->Put(location.x, location.y, location.z, _current);

  // Create a temporary image
  irtkGreyImage tmp(_raster->GetX(), _raster->GetY(), _raster->GetZ());
  tmp = *_raster;

  while (!point_stack.empty()) {
    location = point_stack.top();
    point_stack.pop();
    x = location.x;
    y = location.y;
    z = location.z;
    if (x-1 >= 0) {
      if ((tmp(x-1, y, z) == 0) && (RegionGrowingCriteria(x-1, y, z, lowT, highT))) {
        location.x = x-1;
        location.y = y;
        location.z = z;
        point_stack.push(location);
        tmp(x-1, y, z) = _current;
      }
    }
    if (x+1 < _raster->GetX()) {
      if ((tmp(x+1, y, z) == 0) && (RegionGrowingCriteria(x+1, y, z, lowT, highT))) {
        location.x = x+1;
        location.y = y;
        location.z = z;
        point_stack.push(location);
        tmp(x+1, y, z) = _current;
      }
    }
    if (y-1 >= 0) {
      if ((tmp(x, y-1, z) == 0) && (RegionGrowingCriteria(x, y-1, z, lowT, highT))) {
        location.x = x;
        location.y = y-1;
        location.z = z;
        point_stack.push(location);
        tmp(x, y-1, z) = _current;
      }
    }
    if (y+1 < _raster->GetY()) {
      if ((tmp(x, y+1, z) == 0) && (RegionGrowingCriteria(x, y+1, z, lowT, highT))) {
        location.x = x;
        location.y = y+1;
        location.z = z;
        point_stack.push(location);
        tmp(x, y+1, z) = _current;
      }
    }
  }
  for (z = 0; z < _raster->GetZ(); z++) {
    for (y = 0; y < _raster->GetY(); y++) {
      for (x = 0; x < _raster->GetX(); x++) {
        if ((tmp(x, y, z) == _current) && (_raster->Get(x, y, z) == 0)) {
          _raster->Put(x, y, z, _current);
          _currentSize++;
        }
      }
    }
  }
}

void irtkVoxelContour::RegionGrowing3D(int seedX, int seedY, int seedZ, double lowT, double highT)
{
  int x, y, z;
  irtkLocation location;

  // Create stack
  stack<irtkLocation> point_stack;

  // Push seed location on stack
  location.x = seedX;
  location.y = seedY;
  location.z = seedZ;
  point_stack.push(location);
  _raster->Put(location.x, location.y, location.z, _current);

  // Create a temporary image
  irtkGreyImage tmp(_raster->GetX(), _raster->GetY(), _raster->GetZ());
  tmp = *_raster;

  while (!point_stack.empty()) {
    location = point_stack.top();
    point_stack.pop();
    x = location.x;
    y = location.y;
    z = location.z;
    if (x-1 >= 0) {
      if ((tmp(x-1, y, z) == 0) && (RegionGrowingCriteria(x-1, y, z, lowT, highT))) {
        location.x = x-1;
        location.y = y;
        location.z = z;
        point_stack.push(location);
        tmp(x-1, y, z) = _current;
      }
    }
    if (x+1 < _raster->GetX()) {
      if ((tmp(x+1, y, z) == 0) && (RegionGrowingCriteria(x+1, y, z, lowT, highT))) {
        location.x = x+1;
        location.y = y;
        location.z = z;
        point_stack.push(location);
        tmp(x+1, y, z) = _current;
      }
    }
    if (y-1 >= 0) {
      if ((tmp(x, y-1, z) == 0) && (RegionGrowingCriteria(x, y-1, z, lowT, highT))) {
        location.x = x;
        location.y = y-1;
        location.z = z;
        point_stack.push(location);
        tmp(x, y-1, z) = _current;
      }
    }
    if (y+1 < _raster->GetY()) {
      if ((tmp(x, y+1, z) == 0) && (RegionGrowingCriteria(x, y+1, z, lowT, highT))) {
        location.x = x;
        location.y = y+1;
        location.z = z;
        point_stack.push(location);
        tmp(x, y+1, z) = _current;
      }
    }
    if (z-1 >= 0) {
      if ((tmp(x, y, z-1) == 0) && (RegionGrowingCriteria(x, y, z-1, lowT, highT))) {
        location.x = x;
        location.y = y;
        location.z = z-1;
        point_stack.push(location);
        tmp(x, y, z-1) = _current;
      }
    }
    if (z+1 < _raster->GetZ()) {
      if ((tmp(x, y, z+1) == 0) && (RegionGrowingCriteria(x, y, z+1, lowT, highT))) {
        location.x = x;
        location.y = y;
        location.z = z+1;
        point_stack.push(location);
        tmp(x, y, z+1) = _current;
      }
    }
  }
  for (z = 0; z < _raster->GetZ(); z++) {
    for (y = 0; y < _raster->GetY(); y++) {
      for (x = 0; x < _raster->GetX(); x++) {
        if ((tmp(x, y, z) == _current) && (_raster->Get(x, y, z) == 0)) {
          _raster->Put(x, y, z, _current);
          _currentSize++;
        }
      }
    }
  }
}

