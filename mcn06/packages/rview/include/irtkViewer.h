/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKVIEWER_H

#define _IRTKVIEWER_H

class irtkRView;
class irtkVoxelContour;

class irtkViewer
{

  /// Screen corrdinates
  int _screenX1, _screenY1, _screenX2, _screenY2;

  /// Viewport coordinates
  double _viewportX1, _viewportY1, _viewportX2, _viewportY2;

  /// Pointer to registration viewer
  irtkRView *_rview;

  /// Viewer mode
  irtkViewerMode _viewerMode;

public:

  /// Constructor
  irtkViewer(irtkRView *, irtkViewerMode);

  /// Destructor
  virtual ~irtkViewer();

  /// Draw image viewer
  virtual void DrawImage(irtkColor *);

  /// Draw isolines in image viewer
  virtual void DrawIsolines(irtkGreyImage *, int);

  /// Draw segmentation contours in image viewer
  virtual void DrawSegmentationContour(irtkGreyImage *);

  /// Draw cursor in image viewer
  virtual void DrawCursor(irtkCursorMode mode);

  /// Draw control points
  virtual void DrawPoints();

  /// Draw control points as deformation arrows
  virtual void DrawArrows ();

  /// Draw control points as deformation grid
  virtual void DrawGrid();

  /// Draw landmarks
  void DrawLandmarks(irtkPointSet &, irtkGreyImage *, int = True);

  /// Draw ROI
  void DrawROI(irtkGreyImage *image, double, double, double, double,
               double, double);

#ifdef HAS_VTK
  /// Draw multiple objects
  void DrawObject(vtkPointSet **, irtkGreyImage *, int = False, int = False);

  /// Draw object
  void DrawObject(vtkPointSet *, irtkGreyImage *, int = False, int = False);
#endif

  /// Draw information about L/R, A/P, S/I on the viewer
  void DrawInfo(irtkDisplayMode);

  /// Update using control points
  Bool Update1(irtkGreyImage *, irtkTransformation *);

  /// Update using display resolution
  Bool Update2(irtkGreyImage *, irtkTransformation *);

  /// Update
  Bool Update(irtkGreyImage *, irtkTransformation *);

  /// Get width of viewer
  int GetWidth();

  /// Get height of viewer
  int GetHeight();

  /// Set screen (in normalized coordinates)
  void SetViewport(double, double, double, double);

  /// Get view port (in normalized coordinates)
  void GetViewport(double &, double &, double &, double &);

  /// Set screen coordinates (in pixel coordinates)
  void SetScreen(int, int);

  /// Get view port (in pixel coordinates)
  void GetScreen(int &, int &, int &, int &);

  /// Set viewer mode
  void SetViewerMode(irtkViewerMode);

  /// Get viewer mode
  irtkViewerMode GetViewerMode();

  /// Clipping of a drawable to the viewport
  void  Clip();

};

inline int irtkViewer::GetWidth()
{
  return _screenX2 - _screenX1 + 1;
}

inline int irtkViewer::GetHeight()
{
  return _screenY2 - _screenY1 + 1;
}

inline void irtkViewer::SetViewport(double x1, double y1, double x2, double y2)
{
  _viewportX1 = x1;
  _viewportY1 = y1;
  _viewportX2 = x2;
  _viewportY2 = y2;
}

inline void irtkViewer::GetViewport(double &x1, double &y1, double &x2, double &y2)
{
  x1 = _viewportX1;
  y1 = _viewportY1;
  x2 = _viewportX2;
  y2 = _viewportY2;
}

inline void irtkViewer::SetScreen(int x, int y)
{
  _screenX1 = round(x * _viewportX1);
  _screenY1 = round(y * _viewportY1);
  _screenX2 = round(x * _viewportX2);
  _screenY2 = round(y * _viewportY2);
}

inline void irtkViewer::GetScreen(int &x1, int &y1, int &x2, int &y2)
{
  x1 = _screenX1;
  y1 = _screenY1;
  x2 = _screenX2;
  y2 = _screenY2;
}

inline void irtkViewer::SetViewerMode(irtkViewerMode viewerMode)
{
  _viewerMode = viewerMode;
}

inline irtkViewerMode irtkViewer::GetViewerMode()
{
  return _viewerMode;
}

inline void irtkViewer::Clip()
{
  glViewport(_screenX1, _screenY1, this->GetWidth(), this->GetHeight());
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(_screenX1, _screenX2, _screenY1, _screenY2);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

#endif



