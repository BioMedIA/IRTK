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

#include <irtkRView.h>

// Define the maximum number of control points along each axis we can
// handle here.

#define MaxNumberOfCP 1000

// Define some global variables. These variables should NOT be used by
// any other classes. Perhaps they should be static members.

int     _NumberOfX;
int     _NumberOfY;
double  _AfterX  [MaxNumberOfCP][MaxNumberOfCP];
double  _AfterY  [MaxNumberOfCP][MaxNumberOfCP];
double  _AfterZ  [MaxNumberOfCP][MaxNumberOfCP];
double  _BeforeX [MaxNumberOfCP][MaxNumberOfCP];
double  _BeforeY [MaxNumberOfCP][MaxNumberOfCP];
double  _BeforeZ [MaxNumberOfCP][MaxNumberOfCP];
_Status _CPStatus[MaxNumberOfCP][MaxNumberOfCP];

// Define the default color scheme
#define COLOR_GRID             glColor3f(1, 1, 0)
#define COLOR_ARROWS           glColor3f(1, 1, 0)
#define COLOR_ISOLINES         glColor3f(1, 1, 0)
#define COLOR_CONTOUR          glColor4f(0, 1, 0, 0.5)
#define COLOR_CURSOR           glColor3f(0, 1, 0)
#define COLOR_POINTS_ACTIVE    glColor3f(0, 1, 0)
#define COLOR_POINTS_PASSIVE   glColor3f(0, 0, 1)
#define COLOR_POINTS_UNKNOWN   glColor3f(1, 1, 0)
#define COLOR_TARGET_LANDMARKS glColor3f(1, 0, 0)
#define COLOR_SOURCE_LANDMARKS glColor3f(0, 0, 1)

#ifdef HAS_VTK

// vtk includes
#include <vtkPlane.h>
#include <vtkPolyDataSource.h>
#include <vtkCutter.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkGeometryFilter.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>

// object colour defines
#define COLOR_OBJECT           glColor3f(1, 0, 0)
#endif

GLubyte space[] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};

GLubyte letters[][13] = {
  {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c, 0x18},
  {0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
  {0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
  {0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce, 0xfc},
  {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0, 0xff},
  {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xff},
  {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e},
  {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
  {0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e},
  {0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06},
  {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6, 0xc3},
  {0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0},
  {0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7, 0xc3},
  {0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3, 0xe3},
  {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7, 0x7e},
  {0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
  {0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66, 0x3c},
  {0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe},
  {0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7, 0x7e},
  {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0xff},
  {0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
  {0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
  {0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3},
  {0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
  {0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3},
  {0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03, 0xff}
};

GLuint fontOffset;

// Little helper(s)
void status_glColor(int status)
{
  switch (status) {
  case _Active:
    COLOR_POINTS_ACTIVE;
    break;
  case _Passive:
    COLOR_POINTS_PASSIVE;
    break;
  case _Unknown:
    COLOR_POINTS_UNKNOWN;
    break;
  default:
    break;
  }
}

irtkViewer::irtkViewer(irtkRView *rview, irtkViewerMode viewerMode)
{
  _screenX1   = 0;
  _screenY1   = 0;
  _screenX2   = 0;
  _screenY2   = 0;
  _viewportX1 = 0;
  _viewportY1 = 0;
  _viewportX2 = 0;
  _viewportY2 = 0;

  // Pointer to registration viewer
  _rview = rview;

  // Mode of image viewer
  _viewerMode = viewerMode;
}

irtkViewer::~irtkViewer()
{}

Bool irtkViewer::Update1(irtkGreyImage *image, irtkTransformation *transformation)
{
  double x1, y1, z1, x2, y2, z2, t;
  int i, j, k, i1, j1, k1, i2, j2, k2, index, m, n;

  // Check transformation
  irtkMultiLevelFreeFormTransformation *mffd = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transformation);
  if (mffd == 0) {
    _NumberOfX = 0;
    _NumberOfY = 0;
    return False;
  }

  // Check adaptive FFD
  irtkFreeFormTransformation *affd = (irtkFreeFormTransformation *)mffd->GetLocalTransformation(mffd->NumberOfLevels()-1);
  if ((affd->GetX() >= MaxNumberOfCP) || (affd->GetY() >= MaxNumberOfCP) || (affd->GetZ() >= MaxNumberOfCP)) {
    cerr << "Too many control points" << endl;
    exit(1);
  }

  // Find out time
  t = image->ImageToTime(0);

  // Find out first corner of ROI
  x1 = 0;
  y1 = 0;
  z1 = 0;
  image->ImageToWorld(x1, y1, z1);
  affd->WorldToLattice(x1, y1, z1);
  i1 = round(x1);
  j1 = round(y1);
  k1 = round(z1);

  if (i1 < 0) i1 = 0;
  if (i1 > affd->GetX()-1) i1 = affd->GetX()-1;
  if (j1 < 0) j1 = 0;
  if (j1 > affd->GetY()-1) j1 = affd->GetY()-1;
  if (k1 < 0) k1 = 0;
  if (k1 > affd->GetZ()-1) k1 = affd->GetZ()-1;

  // Find out second corner of ROI
  x2 = image->GetX()-1;
  y2 = image->GetY()-1;
  z2 = 0;
  image->ImageToWorld(x2, y2, z2);
  affd->WorldToLattice(x2, y2, z2);
  i2 = round(x2);
  j2 = round(y2);
  k2 = round(z2);

  if (i2 < 0) i2 = 0;
  if (i2 > affd->GetX()-1) i2 = affd->GetX()-1;
  if (j2 < 0) j2 = 0;
  if (j2 > affd->GetY()-1) j2 = affd->GetY()-1;
  if (k2 < 0) k2 = 0;
  if (k2 > affd->GetZ()-1) k2 = affd->GetZ()-1;

  // Swap if necessary
  if (i1 > i2) swap(i1, i2);
  if (j1 > j2) swap(j1, j2);
  if (k1 > k2) swap(k1, k2);

  switch (_viewerMode) {
  case Viewer_XY:
    _NumberOfX = i2-i1+1;
    _NumberOfY = j2-j1+1;
    break;
  case Viewer_XZ:
    _NumberOfX = i2-i1+1;
    _NumberOfY = k2-k1+1;
    break;
  case Viewer_YZ:
    _NumberOfX = j2-j1+1;
    _NumberOfY = k2-k1+1;
    break;
  }

  for (k = k1; k <= k2; k++) {
    for (j = j1; j <= j2; j++) {
      for (i = i1; i <= i2; i++) {
        // Calculate control points before and after deformation
        switch (_viewerMode) {
        case Viewer_XY:
          m = i-i1;
          n = j-j1;
          break;
        case Viewer_XZ:
          m = i-i1;
          n = k-k1;
          break;
        case Viewer_YZ:
          m = j-j1;
          n = k-k1;
          break;
        default:
          m = i;
          n = j;
          break;
        }

        _BeforeX[m][n] = i;
        _BeforeY[m][n] = j;
        _BeforeZ[m][n] = k;
        index = i * affd->GetY() * affd->GetZ() + j * affd->GetZ() + k;
        affd->LatticeToWorld(_BeforeX[m][n], _BeforeY[m][n], _BeforeZ[m][n]);
        _AfterX[m][n] = _BeforeX[m][n];
        _AfterY[m][n] = _BeforeY[m][n];
        _AfterZ[m][n] = _BeforeZ[m][n];
        mffd->Transform(_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n], t);
        mffd->irtkAffineTransformation::Inverse(_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n]);
        image->WorldToImage(_BeforeX[m][n], _BeforeY[m][n], _BeforeZ[m][n]);
        image->WorldToImage(_AfterX[m][n], _AfterY[m][n], _AfterZ[m][n]);
        if (affd->GetStatus(index) == _Active) {
          _CPStatus[m][n] = _Active;
        } else {
          _CPStatus[m][n] = _Passive;
        }
#ifndef IMPERIAL
        _CPLabel[m][n] = affd->GetLabel(index);
        if (_CPLabel[m][n] > _CPMaxLabel) _CPMaxLabel = _CPLabel[m][n];
        if (_CPLabel[m][n] < _CPMinLabel) _CPMinLabel = _CPLabel[m][n];
#endif

      }
    }
  }
  return True;
}

Bool irtkViewer::Update2(irtkGreyImage *image, irtkTransformation *transformation)
{
  double dx, dy, t;
  int i, j, imax, jmax;

  irtkMultiLevelFreeFormTransformation *mffd = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(transformation);
  if (mffd == NULL) return False;

  dx = _rview->_DisplayDeformationGridResolution;
  dy = _rview->_DisplayDeformationGridResolution;
  _NumberOfX = round((this->GetWidth()  - 40) / (double)dx);
  _NumberOfY = round((this->GetHeight() - 40) / (double)dy);
  imax = _NumberOfX;
  jmax = _NumberOfY;
  dx = (this->GetWidth()  - 40) / (double)_NumberOfX;
  dy = (this->GetHeight() - 40) / (double)_NumberOfY;

  // Find out time
  t = image->ImageToTime(0);

  for (j = 0; j < jmax; j++) {
    for (i = 0; i < imax; i++) {
      _BeforeX[i][j] = i * dx + 20 + dx / 2.0;
      _BeforeY[i][j] = j * dy + 20 + dy / 2.0;
      _BeforeZ[i][j] = 0;
      _AfterX[i][j]  = _BeforeX[i][j];
      _AfterY[i][j]  = _BeforeY[i][j];
      _AfterZ[i][j]  = 0;
      image->ImageToWorld(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j]);
      if (_rview->_sourceTransformInvert == True) {
        mffd->Inverse(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], t);
        mffd->irtkAffineTransformation::Transform(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], t);
      } else {
        mffd->Transform(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], t);
        mffd->irtkAffineTransformation::Inverse(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j], t);
      }
      image->WorldToImage(_AfterX[i][j], _AfterY[i][j], _AfterZ[i][j]);
      _CPStatus[i][j] = _Unknown;
    }
  }
  return True;
}

Bool irtkViewer::Update(irtkGreyImage *image, irtkTransformation *transformation)
{
  if (_rview->_DisplayDeformationGridResolution == 0) {
    return this->Update1(image, transformation);
  } else {
    return this->Update2(image, transformation);
  }
}

void irtkViewer::DrawCursor(irtkCursorMode mode)
{
  int x, y;

  // Set color
  COLOR_CURSOR;

  // calculate width and height
  x = this->GetWidth();
  y = this->GetHeight();

  switch (mode) {
  case CrossHair:
    // Draw cross hair
    glBegin(GL_LINES);
    glVertex2f(_screenX1+x/2-10, _screenY1+y/2);
    glVertex2f(_screenX1+x/2+10, _screenY1+y/2);
    glVertex2f(_screenX1+x/2,    _screenY1+y/2-10);
    glVertex2f(_screenX1+x/2,    _screenY1+y/2+10);
    glEnd();
    break;
  case CursorX:
    // Draw cursor as broken 'X' (+)
    glBegin(GL_LINES);
    glVertex2f(_screenX1+x/2-10, _screenY1+y/2);
    glVertex2f(_screenX1+x/2-3,  _screenY1+y/2);
    glVertex2f(_screenX1+x/2+3,  _screenY1+y/2);
    glVertex2f(_screenX1+x/2+10, _screenY1+y/2);
    glVertex2f(_screenX1+x/2,    _screenY1+y/2-10);
    glVertex2f(_screenX1+x/2,    _screenY1+y/2-3);
    glVertex2f(_screenX1+x/2,    _screenY1+y/2+3);
    glVertex2f(_screenX1+x/2,    _screenY1+y/2+10);
    glEnd();
    break;
  case CursorV:
    // Draw cursor as 'V'
    glBegin(GL_LINES);
    glVertex2f(_screenX1+x/2-4, _screenY1+y/2+10);
    glVertex2f(_screenX1+x/2,    _screenY1+y/2);
    glVertex2f(_screenX1+x/2,    _screenY1+y/2);
    glVertex2f(_screenX1+x/2+4, _screenY1+y/2+10);
    glEnd();
    break;
  case CursorBar:
    // Draw cursor as bar with scales
    break;
  }
}

void irtkViewer::DrawIsolines(irtkGreyImage *image, int value)
{
  int i, j;

  // Set color
  COLOR_ISOLINES;

  glBegin(GL_LINES);

  for (j = 0; j < this->GetHeight()-1; j++) {
    for (i = 0; i < this->GetWidth()-1; i++) {
      if (((image->Get(i, j, 0)   <= value) &&
           (image->Get(i+1, j, 0)  > value)) ||
          ((image->Get(i, j, 0)    > value) &&
           (image->Get(i+1, j, 0) <= value))) {
        glVertex2f(_screenX1+i+0.5, _screenY1+j-0.5);
        glVertex2f(_screenX1+i+0.5, _screenY1+j+0.5);
      }
      if (((image->Get(i, j, 0)   <= value) &&
           (image->Get(i, j+1, 0)  > value)) ||
          ((image->Get(i, j, 0)    > value) &&
           (image->Get(i, j+1, 0) <= value))) {
        glVertex2f(_screenX1+i+0.5, _screenY1+j+0.5);
        glVertex2f(_screenX1+i-0.5, _screenY1+j+0.5);
      }
    }
  }
  glEnd();
}

void irtkViewer::DrawSegmentationContour(irtkGreyImage *image)
{
  int i, j;
  unsigned char r, g, b;

  for (j = 1; j < this->GetHeight()-1; j++) {
    for (i = 1; i < this->GetWidth()-1; i++) {
      if ((image->Get(i, j, 0) > 0) && (_rview->_segmentTable->_entry[image->Get(i, j, 0)]._visible == True)) {
        r = _rview->_segmentTable->_entry[image->Get(i, j, 0)]._color.r;
        g = _rview->_segmentTable->_entry[image->Get(i, j, 0)]._color.g;
        b = _rview->_segmentTable->_entry[image->Get(i, j, 0)]._color.b;
        if (image->Get(i, j, 0) != image->Get(i+1, j, 0)) {
          glColor3ub(r, g, b);
          glBegin(GL_LINES);
          glVertex2f(_screenX1+i, _screenY1+j-0.5);
          glVertex2f(_screenX1+i, _screenY1+j+0.5);
          glEnd();
        }
        if (image->Get(i, j, 0) != image->Get(i-1, j, 0)) {
          glColor3ub(r, g, b);
          glBegin(GL_LINES);
          glVertex2f(_screenX1+i, _screenY1+j-0.5);
          glVertex2f(_screenX1+i, _screenY1+j+0.5);
          glEnd();
        }
        if (image->Get(i, j, 0) != image->Get(i, j+1, 0)) {
          glColor3ub(r, g, b);
          glBegin(GL_LINES);
          glVertex2f(_screenX1+i+0.5, _screenY1+j);
          glVertex2f(_screenX1+i-0.5, _screenY1+j);
          glEnd();
        }
        if (image->Get(i, j, 0) != image->Get(i, j-1, 0)) {
          glColor3ub(r, g, b);
          glBegin(GL_LINES);
          glVertex2f(_screenX1+i+0.5, _screenY1+j);
          glVertex2f(_screenX1+i-0.5, _screenY1+j);
          glEnd();
        }
      }
    }
  }
}

void irtkViewer::DrawGrid()
{
  int i, j;

  // Set color
  COLOR_GRID;

  glBegin(GL_LINES);
  for (j = 0; j < _NumberOfY; j++) {
    for (i = 0; i < _NumberOfX-1; i++) {
      glVertex2f(_screenX1+_AfterX[i][j],   _screenY1+_AfterY[i][j]);
      glVertex2f(_screenX1+_AfterX[i+1][j], _screenY1+_AfterY[i+1][j]);
    }
  }
  for (j = 0; j < _NumberOfY-1; j++) {
    for (i = 0; i < _NumberOfX; i++) {
      glVertex2f(_screenX1+_AfterX[i][j],   _screenY1+_AfterY[i][j]);
      glVertex2f(_screenX1+_AfterX[i][j+1], _screenY1+_AfterY[i][j+1]);
    }
  }
  glEnd();
}

void irtkViewer::DrawArrows()
{
  int i, j;

  // Set color
  COLOR_ARROWS;

  for (j = 0; j < _NumberOfY; j++) {
    for (i = 0; i < _NumberOfX; i++) {
      glBegin(GL_LINES);
      glVertex2f(_screenX1+_BeforeX[i][j], _screenY1+_BeforeY[i][j]);
      glVertex2f(_screenX1+_AfterX[i][j],  _screenY1+_AfterY[i][j]);
      glEnd();
    }
  }
}

void irtkViewer::DrawPoints()
{
  int i, j;

  // Adjust pointsize
  glPointSize((GLfloat)3);

  // Draw active and passive points
  glBegin(GL_POINTS);
  for (j = 0; j < _NumberOfY; j++) {
    for (i = 0; i < _NumberOfX; i++) {
      // Set color
      status_glColor(_CPStatus[i][j]);
      // Draw point
      glVertex2f(_screenX1+_BeforeX[i][j], _screenY1+_BeforeY[i][j]);
    }
  }
  glEnd();
}

void irtkViewer::DrawLandmarks(irtkPointSet &landmarks, irtkGreyImage *image,
                               int bTarget)
{
  int i;
  irtkPoint p;

  // Adjust pointsize and colour
  if (bTarget == True) {
    COLOR_TARGET_LANDMARKS;
  } else {
    COLOR_SOURCE_LANDMARKS;
  }

  // Draw landmarks
  for (i = 0; i < landmarks.Size(); i++) {

    // Get point
    p = landmarks(i);

    if (bTarget == False) {
      // Transform point
      _rview->_sourceTransform->Inverse(p._x, p._y, p._z);
    }

    // Draw point
    if (image->IsInFOV(p._x, p._y, p._z) == True) {
      image->WorldToImage(p);
      glBegin(GL_LINES);
      glVertex2f(_screenX1+p._x-8, _screenY1+p._y);
      glVertex2f(_screenX1+p._x+8, _screenY1+p._y);
      glVertex2f(_screenX1+p._x,   _screenY1+p._y-8);
      glVertex2f(_screenX1+p._x,   _screenY1+p._y+8);
      glEnd();
    }
  }
}

void irtkViewer::DrawImage(irtkColor *drawable)
{
  // Set raster position
  glRasterPos2f(_screenX1, _screenY1);

  // Draw pixelmap
  glDrawPixels(this->GetWidth(), this->GetHeight(), GL_RGB, GL_UNSIGNED_BYTE,
               drawable);
}

void irtkViewer::DrawROI(irtkGreyImage *image,
                         double x1, double y1, double z1,
                         double x2, double y2, double z2)
{
  image->WorldToImage(x1, y1, z1);
  image->WorldToImage(x2, y2, z2);
  glColor4f(1, 0, 0, 1);
  glBegin(GL_POLYGON);
  glVertex2f(_screenX1+x1-2, _screenY1+y1-2);
  glVertex2f(_screenX1+x1+2, _screenY1+y1-2);
  glVertex2f(_screenX1+x1+2, _screenY1+y1+2);
  glVertex2f(_screenX1+x1-2, _screenY1+y1+2);
  glEnd();
  glColor4f(0, 1, 0, 1);
  glBegin(GL_POLYGON);
  glVertex2f(_screenX1+x2-2, _screenY1+y2-2);
  glVertex2f(_screenX1+x2+2, _screenY1+y2-2);
  glVertex2f(_screenX1+x2+2, _screenY1+y2+2);
  glVertex2f(_screenX1+x2-2, _screenY1+y2+2);
  glEnd();
  glColor4f(1, 1, 0, 0.5);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glBegin(GL_POLYGON);
  glVertex2f(_screenX1+x1, _screenY1+y1);
  glVertex2f(_screenX1+x1, _screenY1+y2);
  glVertex2f(_screenX1+x2, _screenY1+y2);
  glVertex2f(_screenX1+x2, _screenY1+y1);
  glEnd();
  glDisable(GL_BLEND);
}


#ifdef HAS_VTK

void irtkViewer::DrawObject(vtkPointSet **object, irtkGreyImage *image,
         int _DisplayObjectWarp,
         int _DisplayObjectGrid)
{
  int i;

  for (i = 0; i < _rview->_NoOfObjects; i++) {
    glColor3ub(_rview->GetObjectLookupTable()->lookupTable[i].r,
           _rview->GetObjectLookupTable()->lookupTable[i].g,
           _rview->GetObjectLookupTable()->lookupTable[i].b);

    this->DrawObject(object[i], image, _DisplayObjectWarp, _DisplayObjectGrid);
  }
}

void irtkViewer::DrawObject(vtkPointSet *points, irtkGreyImage *image, int _DisplayObjectWarp, int _DisplayObjectGrid)
{
  int i, j;
  double point[3], bounds[6], origin[3], normal[3];
  irtkPoint p;
  static vtkPlane *plane   = vtkPlane::New();
  static vtkCutter *cutter = vtkCutter::New();

  if (points != NULL) {

    points->GetCenter(origin);
    points->GetBounds(bounds);
    switch (_viewerMode) {
    case Viewer_XY:
      origin[2] = (int)round(image->GetOrigin()._z);
      normal[0] = 0;
      normal[1] = 0;
      normal[2] = (-origin[1]*(bounds[0]-origin[0]) +
                   origin[0]*(bounds[1]-origin[1]));
      break;
    case Viewer_YZ:
      origin[0] = (int)round(image->GetOrigin()._x);
      normal[0] = (-origin[2]*(bounds[1]-origin[1]) +
                   origin[1]*(bounds[2]-origin[2]));
      normal[1] = 0;
      normal[2] = 0;
      break;
    case Viewer_XZ:
      origin[1] = (int)round(image->GetOrigin()._y);
      normal[0] = 0;
      normal[1] = (-origin[2]*(bounds[0]-origin[0]) +
                   origin[0]*(bounds[2]-origin[2]));
      normal[2] = 0;
      break;
    }

    // Set up plane and cutter
    plane->SetOrigin(origin);
    plane->SetNormal(normal);
    cutter->SetCutFunction(plane);
    cutter->SetInput(points);

    // Reslice object
    cutter->Modified();
    cutter->Update();

    // Loop over cells
    for (i = 0; i < cutter->GetOutput()->GetNumberOfCells(); i++) {
      // Get pointIds from cell
      vtkIdList *ids = cutter->GetOutput()->GetCell(i)->GetPointIds();
      irtkPointSet pset;
      for (j = 0; j < ids->GetNumberOfIds(); j++) {

        // Get point from cell
        cutter->GetOutput()->GetPoints()->GetPoint(ids->GetId(j), point);

        // Get point in image coordinates
        image->WorldToImage(point[0], point[1], point[2]);
        switch (_viewerMode) {
        case Viewer_XY:
          p = irtkPoint(point[0], point[1], 0);
          break;
        case Viewer_YZ:
          p = irtkPoint(point[0], point[1], 0);
          break;
        case Viewer_XZ:
          p = irtkPoint(point[0], point[1], 0);
          break;
        }
        pset.Add(p);
      }

      // Now draw
      for (j = 0; j < pset.Size(); j++) {
        glBegin(GL_LINES);
        glVertex2f(_screenX1+pset(j)._x,
                   _screenY1+pset(j)._y);
        glVertex2f(_screenX1+pset((j+1)%pset.Size())._x,
                   _screenY1+pset((j+1)%pset.Size())._y);
        glEnd();
      }
    }
  }
}

#endif

void irtkViewer::DrawInfo(irtkDisplayMode m)
{
  int x, y;
  static int first = True;

  if (first == True) {
    GLuint i, j;
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    fontOffset = glGenLists (128);
    for (i = 0,j = 'A'; i < 26; i++,j++) {
      glNewList(fontOffset + j, GL_COMPILE);
      glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, letters[i]);
      glEndList();
    }
    glNewList(fontOffset + ' ', GL_COMPILE);
    glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, space);
    glEndList();
    first = False;
  }

  if (m == Native) return;

  // calculate width and height
  x = this->GetWidth();
  y = this->GetHeight();

  switch (_viewerMode) {
  case Viewer_XY:

    if (m == Neurological) {

      // Draw axis labels
      glRasterPos2i(_screenX1+5, _screenY1+y/2-5);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(strlen("L"), GL_UNSIGNED_BYTE, (GLubyte *) "R");
      glPopAttrib ();

      glRasterPos2i(_screenX1+x-15, _screenY1+y/2-5);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(strlen("R"), GL_UNSIGNED_BYTE, (GLubyte *) "L");
      glPopAttrib ();

    } else {

      // Draw axis labels
      glRasterPos2i(_screenX1+5, _screenY1+y/2-5);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(strlen("L"), GL_UNSIGNED_BYTE, (GLubyte *) "L");
      glPopAttrib ();

      glRasterPos2i(_screenX1+x-15, _screenY1+y/2-5);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(strlen("R"), GL_UNSIGNED_BYTE, (GLubyte *) "R");
      glPopAttrib ();

    }

    glRasterPos2i(_screenX1+x/2-5, _screenY1+5);
    glPushAttrib (GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists(strlen("P"), GL_UNSIGNED_BYTE, (GLubyte *) "P");
    glPopAttrib ();

    glRasterPos2i(_screenX1+x/2-5, _screenY1+y-15);
    glPushAttrib (GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists(strlen("A"), GL_UNSIGNED_BYTE, (GLubyte *) "A");
    glPopAttrib ();
    break;

  case Viewer_XZ:

    if (m == Neurological) {

      // Draw axis labels
      glRasterPos2i(_screenX1+5, _screenY1+y/2-5);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(strlen("R"), GL_UNSIGNED_BYTE, (GLubyte *) "R");
      glPopAttrib ();

      glRasterPos2i(_screenX1+x-15, _screenY1+y/2-5);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(strlen("L"), GL_UNSIGNED_BYTE, (GLubyte *) "L");
      glPopAttrib ();

    } else {

      // Draw axis labels
      glRasterPos2i(_screenX1+5, _screenY1+y/2-5);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(strlen("R"), GL_UNSIGNED_BYTE, (GLubyte *) "L");
      glPopAttrib ();

      glRasterPos2i(_screenX1+x-15, _screenY1+y/2-5);
      glPushAttrib (GL_LIST_BIT);
      glListBase(fontOffset);
      glCallLists(strlen("L"), GL_UNSIGNED_BYTE, (GLubyte *) "R");
      glPopAttrib ();

    }

    glRasterPos2i(_screenX1+x/2-5, _screenY1+5);
    glPushAttrib (GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists(strlen("I"), GL_UNSIGNED_BYTE, (GLubyte *) "I");
    glPopAttrib ();

    glRasterPos2i(_screenX1+x/2-5, _screenY1+y-15);
    glPushAttrib (GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists(strlen("S"), GL_UNSIGNED_BYTE, (GLubyte *) "S");
    glPopAttrib ();
    break;

  case Viewer_YZ:

    // Draw axis labels
    glRasterPos2i(_screenX1+5, _screenY1+y/2-5);
    glPushAttrib (GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists(strlen("P"), GL_UNSIGNED_BYTE, (GLubyte *) "P");
    glPopAttrib ();

    glRasterPos2i(_screenX1+x-15, _screenY1+y/2-5);
    glPushAttrib (GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists(strlen("A"), GL_UNSIGNED_BYTE, (GLubyte *) "A");
    glPopAttrib ();

    glRasterPos2i(_screenX1+x/2-5, _screenY1+5);
    glPushAttrib (GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists(strlen("I"), GL_UNSIGNED_BYTE, (GLubyte *) "I");
    glPopAttrib ();

    glRasterPos2i(_screenX1+x/2-5, _screenY1+y-15);
    glPushAttrib (GL_LIST_BIT);
    glListBase(fontOffset);
    glCallLists(strlen("S"), GL_UNSIGNED_BYTE, (GLubyte *) "S");
    glPopAttrib ();
    break;
  }
}

