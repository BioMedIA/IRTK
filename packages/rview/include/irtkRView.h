/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKRVIEW_H

#define _IRTKRVIEW_H

typedef enum { View_A,
               View_B,
               View_Checkerboard,
               View_Subtraction,
               View_HShutter,
               View_VShutter,
               View_AoverB,
               View_BoverA
             } irtkRViewMode;

typedef enum { NoneDef,
               Displacement,
               Jacobian,
               Jacobian_Expansion,
               Jacobian_Contraction
             } irtkDeformationProperty;

typedef enum { Viewer_XY, Viewer_XZ, Viewer_YZ, Viewer_None } irtkViewerMode;

typedef enum { CrossHair, CursorX, CursorV, CursorBar } irtkCursorMode;

typedef enum { Neurological, Radiological, Native } irtkDisplayMode;

typedef enum { RegionGrowing2D, RegionGrowing3D } irtkRegionGrowingMode;

#define DEFORMATION_DISPLACEMENT_MIN 0
#define DEFORMATION_DISPLACEMENT_MAX 100
#define DEFORMATION_JACOBIAN_MIN 10
#define DEFORMATION_JACOBIAN_MAX 400

// Function keys
#define KEY_F1                     1
#define KEY_F2                     2
#define KEY_F3                     3
#define KEY_F4                     4
#define KEY_F5                     5
#define KEY_F6                     6
#define KEY_F7                     7
#define KEY_F8                     8
#define KEY_F9                     9
#define KEY_F10                    10
#define KEY_F11                    11
#define KEY_F12                    12

// Directional keys
#define KEY_LEFT                   100
#define KEY_UP                     101
#define KEY_RIGHT                  102
#define KEY_DOWN                   103
#define KEY_PAGE_UP                104

// Mouse buttons.
#define MOUSE_LEFT_BUTTON                0
#define MOUSE_MIDDLE_BUTTON              1
#define MOUSE_RIGHT_BUTTON               2

// Mouse button  state.
#define MOUSE_DOWN                       0
#define MOUSE_UP                         1

#ifndef IMPERIAL

#include <irtkTransformationCollection.h>

#define AFFDTransformation "irtkAdaptiveFreeFormTransformation"
#define MFFDTransformation "irtkTransformationCollection"

#define irtkBSplineFreeFormTransformation irtkAdaptiveFreeFormTransformation
#define irtkMFreeFormTransformation irtkTransformationCollection
#define irtkMFreeFormRegistration   irtkAdaptiveFreeFormRegistration
#define irtkMFreeFormRegistration2D irtkAdaptiveFreeFormRegistration2D

#endif

#include <list>

#ifdef HAS_VTK

#include <vtkPolyData.h>
#include <vtkPointSet.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#endif

#define MAX_SEGMENTS 256

#define MAX_NUMBER_OF_OBJECTS 40

#include <irtkSegmentTable.h>

#include <irtkLookupTable.h>
#include <irtkViewer.h>
#include <irtkRViewConfig.h>
#include <irtkHistogramWindow.h>
#include <irtkVoxelContour.h>

class irtkVoxelContour;

class irtkRView
{

protected:

  /// Friends
  friend class irtkViewer;
  friend class irtkLookupTable;
  friend class irtkVoxelContour;
  friend class irtkSegmentationEditor;

  /// Number of image viewers
  int _NoOfViewers;

  /// Image viewer for target image
  irtkViewer **_viewer;

  /// Target image
  irtkImage *_targetImage;

  /// Source image
  irtkImage *_sourceImage;

  /// Segmentation image
  irtkGreyImage *_segmentationImage;

  /// Segment Table
  irtkSegmentTable *_segmentTable;

  /// Color Table
  irtkColor *_segmentColorTable;

  /// Transformation for reslicing of target image
  irtkTransformation *_targetTransform;

  /// Transformation for reslicing of source image
  irtkTransformation *_sourceTransform;

  /// Transformation for reslicing of segmentation image
  irtkTransformation *_segmentationTransform;

  /// Transformation for reslicing of selection image
  irtkTransformation *_selectionTransform;

  /// Transformation filter for reslicing of target image
  irtkImageTransformation **_targetTransformFilter;

  /// Transformation filter for reslicing of source image
  irtkImageTransformation **_sourceTransformFilter;

  /// Transformation filter for reslicing of segmentation image
  irtkImageTransformation **_segmentationTransformFilter;

  /// Transformation filter for reslicing of selection image
  irtkImageTransformation **_selectionTransformFilter;

  /// Target image
  irtkGreyImage **_targetImageOutput;

  /// Source image
  irtkGreyImage **_sourceImageOutput;

  /// Segmentation image
  irtkGreyImage **_segmentationImageOutput;

  /// Selection image
  irtkGreyImage **_selectionImageOutput;

  /// Target landmarks (pointset)
  irtkPointSet _targetLandmarks;

  /// Source landmarks (pointset)
  irtkPointSet _sourceLandmarks;

  /// Contour
  irtkVoxelContour _voxelContour;

  /// Contour viewer
  int _contourViewer;

  /// Contour viewer mode
  irtkViewerMode _contourViewerMode;

#ifdef HAS_VTK
  /// Number of vtk objects
  int _NoOfObjects;

  /// Object (surface or mesh)
  vtkPointSet *_Object[MAX_NUMBER_OF_OBJECTS];

  /// Flag for display object as a movie
  int _ObjectMovie;

#endif

  /// Combined source and target images in OpenGL format
  irtkColor **_drawable;

  /// Color lookup table for target image
  irtkLookupTable *_targetLookupTable;

  /// Color lookup table for source image
  irtkLookupTable *_sourceLookupTable;

  /// Color lookup table for subtraction of target and source image
  irtkLookupTable *_subtractionLookupTable;

  /// Target value range
  double _targetMin, _targetMax;

  /// Source value range
  double _sourceMin, _sourceMax;

  /// Subtraction value range
  double _subtractionMin, _subtractionMax;

  /// Target display value range
  double _targetDisplayMin, _targetDisplayMax;
  
  /// Source display value range
  double _sourceDisplayMin, _sourceDisplayMax;

  /// Subtraction display value range
  double _subtractionDisplayMin, _subtractionDisplayMax;

  /// Target frame
  int _targetFrame;

  /// Source frame
  int _sourceFrame;

  /// Flag to indicate whether target image must be updated
  int _targetUpdate;

  /// Flag to indicate whether source image must be updated
  int _sourceUpdate;

  /// Flag to indicate whether source image must be updated
  int _segmentationUpdate;

  /// Flag to indicate whether selection image must be updated
  int _selectionUpdate;

  /// Width of viewer  (in pixels)
  int _screenX;

  /// Height of viewer (in pixels)
  int _screenY;

  /// Display origin (in mm)
  double _origin_x;

  /// Display origin (in mm)
  double _origin_y;

  /// Display origin (in mm)
  double _origin_z;

  /// Display resolution
  double _resolution;

  /// Display axes
  double _xaxis[3], _yaxis[3], _zaxis[3];

  /// Region of interest
  double _x1, _y1, _z1, _x2, _y2, _z2;

  /// Current mouse
  /// Interpolator for target image
  irtkImageFunction *_targetInterpolator;

  /// Interpolator for source image
  irtkImageFunction *_sourceInterpolator;

  /// Interpolator for segmentation image
  irtkImageFunction *_segmentationInterpolator;

  /// Interpolator for selection image
  irtkImageFunction *_selectionInterpolator;

  /// Flag whether transformation for reslicing of source image should be applied
  bool _sourceTransformApply;

  /// Flag whether transformation for reslicing of source image should be inverted
  bool _sourceTransformInvert;

  /// Display viewing mix in shutter viewing mode
  double _viewMix;

  /// Flag for rview mode
  irtkRViewMode _viewMode;

  /// Flag for configuration mode
  irtkConfigViewerMode _configMode;

  /// Flag for display labels of segmentation image
  int _DisplaySegmentationLabels;

  /// Flag for display contours of segmentation image
  int _DisplaySegmentationContours;

  /// Flag for segmentation mode
  int _SegmentationMode;

  /// Parameter for paintbrush width
  int _PaintBrushWidth;

  /// Minimum region growing threshold
  int _RegionGrowingThresholdMin;

  /// Maximum region growing threshold
  int _RegionGrowingThresholdMax;

  /// Deformation property
  irtkDeformationProperty _DeformationProperty;

  /// Deformation blending
  double _DeformationBlending;

  /// Flag for display orientation
  irtkDisplayMode _DisplayMode;

  /// Flag for line thickness
  double _LineThickness;

  /// Flag for display of isolines from target image
  int _DisplayTargetContour;

  /// Flag for display of isolines from source image
  int _DisplaySourceContour;

  /// Flag for snap to grid
  int _SnapToGrid;

  /// Flag for display of cross hair
  int _DisplayCursor;

  /// Mode for cursor display
  irtkCursorMode _CursorMode;

  /// Flag for display of axis labels
  int _DisplayAxisLabels;

  /// Resolution for display of deformation grid
  int _DisplayDeformationGridResolution;

  /// Flag for display of deformation grid
  int _DisplayDeformationGrid;

  /// Flag for display of deformation points
  int _DisplayDeformationPoints;

  /// Flag for display of deformation arrows
  int _DisplayDeformationArrows;

  /// Flag for display of landmarks
  int _DisplayLandmarks;

  /// Flag for display of ROI
  int _DisplayROI;

  /// Flag for track of tag using gravity window
  int _TrackTAG;

  /// Flag for display of tag grid pattern
  int _ViewTAG;

  /// Flag for flipping X coordinates
  int _FlipX;

  /// Flag for flipping Y coordinates
  int _FlipY;

  /// Flag for flipping Z coordinates
  int _FlipZ;

  /// Current mouse position
  int _mouseX, _mouseY, _mouseZ;

  /// Current intensity at mouse position
  double _mouseTargetIntensity;

  /// Current viewer in which the mouse is
  int _mouseViewer;

  /// Region growing mode
  irtkRegionGrowingMode _regionGrowingMode;

#ifdef HAS_VTK
  /// Flag for display of object
  int _DisplayObject;

  /// Flag for warping of object
  int _DisplayObjectWarp;

  /// Flag for display of object grid
  int _DisplayObjectGrid;
#endif

public:

  /// Constructor
  irtkRView(int, int);

  /// Destructor
  virtual ~irtkRView();

  /// Render
  void Draw();

  /// Render offscreen
  void DrawOffscreen(char *);

  /// Update registration viewer
  void Update();

  /// Set update of source transformation to on
  void SourceUpdateOn();

  /// Set update of segmentation transformation to on
  void SegmentationUpdateOn();

  /// Resize registration viewer
  void Resize(int, int);

  /// Initialize registration viewer
  virtual void Initialize();

  /// Configure registration viewer
  virtual void Configure(irtkRViewConfig []);

  /// Read configuration
  virtual void Read(char *);

  /// Write configuration
  virtual void Write(char *);

  /// Read target image
  virtual void ReadTarget(char *);

  /// Read target image sequence
  virtual void ReadTarget(int, char **);

  /// Read source image
  virtual void ReadSource(char *);

  /// Read source image sequence
  virtual void ReadSource(int, char **);

  /// Read segmentation image
  virtual void ReadSegmentation(char *);

  /// Write target image
  virtual void WriteTarget(char *);

  /// Write source image
  virtual void WriteSource(char *);

  /// Write segmentation image
  virtual void WriteSegmentation(char *);

  /// Read transformation
  virtual void ReadTransformation(char *);

  /// Write transformation
  virtual void WriteTransformation(char *);

  /// Read target landmarks
  virtual void ReadTargetLandmarks(char *);

  /// Read source landmarks
  virtual void ReadSourceLandmarks(char *);

  /// Write target landmarks
  virtual void WriteTargetLandmarks(char *);

  /// Write sourcelandmarks
  virtual void WriteSourceLandmarks(char *);

#ifdef HAS_VTK
  /// Read object
  virtual void ReadObject(const char *);

  /// Remove object
  virtual void RemoveObject();

  /// Return object
  virtual vtkPointSet *GetObject(int);

  /// Turn display of object movie
  void ObjectMovieOn();

  /// Turn display of object movie
  void ObjectMovieOff();

  /// Return display of object movie
  int GetObjectMovie();
#endif

  /// Get width of registration viewer (in pixels)
  int GetWidth();

  /// Get height of registration viewer (in pixels)
  int GetHeight();

  /// Set display resolution (in mm)
  void   SetResolution(double);

  /// Get display resolution (in mm)
  double GetResolution();

  /// Set display origin (in mm)
  void SetOrigin(double,   double,   double);

  /// Set display origin (in display pixels)
  void SetOrigin(int, int);

  /// Get display origin (in mm)
  void GetOrigin(double &, double &, double &);

  /// Set target frame
  void SetTargetFrame(int);

  /// Get target frame
  int GetTargetFrame();

  /// Set source frame
  void SetSourceFrame(int);

  /// Get source frame
  int GetSourceFrame();

  /// Get minimum target intensity
  double GetTargetMin();
  
  /// Get maximum target intensity
  double GetTargetMax();
  
  /// Get minimum source intensity
  double GetSourceMin();
  
  /// Get maximum source intensity
  double GetSourceMax();
  
  /// Get minimum subtraction intensity
  double GetSubtractionMin();
  
  /// Get maximum subtraction intensity
  double GetSubtractionMax();
  
  /// Set ROI to default parameters
  void ResetROI();

  // Update left-upper corner of ROI
  void UpdateROI1(int, int);

  // Update right-lower corner of ROI
  void UpdateROI2(int, int);

  /// Add a point to the contour
  void AddContour(int, int, ContourMode mode);

  /// Undo adding last part contour
  void UndoContour();

  /// Delete current contour
  void ClearContour();

  /// Fill contour
  void FillContour(int, int);

  /// Fills an area defined by contour
  void FillArea(int, int);

  /// Region growing for contour
  void RegionGrowContour(int, int);

  /// Set ROI
  void SetROI(double, double, double, double, double, double);

  /// Get ROI
  void GetROI(double &, double &, double &, double &, double &, double &);

  /// Get an information string
  void GetInfoText(char *, char *, char *, char *, char *);

  /// Current mouse position
  void MousePosition(int, int);

  /// Current mouse wheel
  void MouseWheel(int, int, int);

  /// Set glLine thickness
  void SetLineThickness(double value);

  /// Get glLine thickness
  double GetLineThickness();

  /// Get an information string about the transformation level
  void GetTransformationText(list<char *> &);

  /// Turn iso-contours extracted from target image on
  void DisplayTargetContoursOn();

  /// Turn iso-contours extracted from target image off
  void DisplayTargetContoursOff();

  /// Return display of target iso-contours
  int GetDisplayTargetContours();

  /// Turn iso-contours extracted from source image on
  void DisplaySourceContoursOn();

  /// Turn iso-contours extracted from source image on
  void DisplaySourceContoursOff();

  /// Return display of source iso-contours
  int GetDisplaySourceContours();

  /// Return display mode
  irtkDisplayMode GetDisplayMode();

  /// Set display mode
  void SetDisplayMode(irtkDisplayMode mode);

  /// Turn snap to grid on
  void SnapToGridOn();

  /// Turn snap to grid off
  void SnapToGridOff();

  /// Return snap to grid
  int GetSnapToGrid();

  /// Turn axis labels on
  void DisplayAxisLabelsOn();

  /// Turn axis labels off
  void DisplayAxisLabelsOff();

  /// Return display of cross hair
  int GetDisplayCursor();

  /// Return cursor mode
  irtkCursorMode GetCursorMode();

  /// Set cursor mode
  void SetCursorMode(irtkCursorMode mode);

  /// Return minimum display intensity of target image
  double GetDisplayMinTarget();

  /// Sets minimum display intensity of target image
  void SetDisplayMinTarget(double);

  /// Return maximum display intensity of target image
  double GetDisplayMaxTarget();

  /// Sets maximum display intensity of target image
  void SetDisplayMaxTarget(double);

  /// Return minimum display intensity of source image
  double GetDisplayMinSource();

  /// Sets minimum display intensity of source image
  void SetDisplayMinSource(double);

  /// Return maximum display intensity of source image
  double GetDisplayMaxSource();

  /// Sets maximum display intensity of source image
  void SetDisplayMaxSource(double);

  /// Return minimum display intensity of subtraction image
  double GetDisplayMinSubtraction();

  /// Sets minimum display intensity of subtraction image
  void SetDisplayMinSubtraction(double);

  /// Return maximum display intensity of subtraction image
  double GetDisplayMaxSubtraction();

  /// Sets maximum display intensity of subtraction image
  void SetDisplayMaxSubtraction(double);

  /// Return minimum display intensity of deformation 
  double GetDisplayMinDeformation();

  /// Sets minimum display intensity of deformation 
  void SetDisplayMinDeformation(double);

  /// Return maximum display intensity of deformation
  double GetDisplayMaxDeformation();

  /// Sets maximum display intensity of deformation
  void SetDisplayMaxDeformation(double);

  /// Turn cross hair on
  void DisplayCursorOn();

  /// Turn cross hair off
  void DisplayCursorOff();

  /// Turn display of deformation grid on
  void DisplayDeformationGridOn();

  /// Turn display of deformation grid off
  void DisplayDeformationGridOff();

  /// Return display of deformation grid
  int GetDisplayDeformationGrid();

  /// Return display of deformation grid
  int GetDisplayDeformationGridResolution();

  /// Sets display resolution of deformation grid
  void SetDisplayDeformationGridResolution(int);

  /// Turn display of deformation control points on
  void DisplayDeformationPointsOn();

  /// Turn display of deformation control points off
  void DisplayDeformationPointsOff();

  /// Return display of deformation points
  int GetDisplayDeformationPoints();

  /// Turn display of ROI
  void DisplayROIOn();

  /// Turn display of ROI
  void DisplayROIOff();

  /// Return display of ROI
  int GetDisplayROI();

  /// Turn display of ROI
  void TrackTAGOn();

  /// Turn display of ROI
  void TrackTAGOff();

  /// Return display of ROI
  int GetTrackTAG();

    /// Turn display of ROI
  void ViewTAGOn();

  /// Turn display of ROI
  void ViewTAGOff();

  /// Return display of ROI
  int GetViewTAG();

  /// Turn display of deformation arrows on
  void DisplayDeformationArrowsOn();

  /// Turn display of deformation arrows off
  void DisplayDeformationArrowsOff();

  /// Return display of deformation arrows
  int GetDisplayDeformationArrows();

  /// Return display of relative deformations
  int GetDisplayDeformationRelative();

  /// Turn display of landmarks on
  void DisplayLandmarksOn();

  /// Turn display of landmarks off
  void DisplayLandmarksOff();

  /// Return display of landmarks
  int GetDisplayLandmarks();

#ifdef HAS_VTK
  /// Turn display of object on
  void DisplayObjectOn();

  /// Turn display of object off
  void DisplayObjectOff();

  /// Return display of object
  int GetDisplayObject();

  /// Turn warping of object on
  void DisplayObjectWarpOn();

  /// Turn warping of object off
  void DisplayObjectWarpOff();

  /// Return warping of object
  int GetDisplayObjectWarp();

  /// Turn grid display of object on
  void DisplayObjectGridOn();

  /// Turn grid display of object off
  void DisplayObjectGridOff();

  /// Return grid display of object
  int GetDisplayObjectGrid();
#endif

  /// Set region growing mode
  void SetRegionGrowingMode(irtkRegionGrowingMode);

  /// Get region growing mode
  irtkRegionGrowingMode GetRegionGrowingMode();

  /// Flip X on
  void FlipXOn();

  /// Flip X off
  void FlipXOff();

  /// Flip Y on
  void FlipYOn();

  /// Flip Y off
  void FlipYOff();

  /// Flip Z on
  void FlipZOn();

  /// Flip Z off
  void FlipZOff();

  /// Set the viewing mix for target and source in shutter viewing mode
  void   SetViewMix(double);

  /// Get the viewing mix for target and source in shutter viewing mode
  double GetViewMix();

  /// Get viewing mode for registration viewer
  irtkRViewMode GetViewMode();

  /// Set viewing mode for registration viewer
  void SetViewMode(irtkRViewMode);

  /// Get configuration mode for viewer
  irtkConfigViewerMode GetConfigMode();

  /// Set configuration mode for viewer
  void SetConfigMode(irtkConfigViewerMode );

  /// Set interpolation mode for target image
  void   SetTargetInterpolationMode(irtkInterpolationMode);

  /// Get interpolation mode for target image
  irtkInterpolationMode GetTargetInterpolationMode();

  /// Set interpolation mode for source image
  void   SetSourceInterpolationMode(irtkInterpolationMode);

  /// Get interpolation model fo target image
  irtkInterpolationMode GetSourceInterpolationMode();

  /// Set transformation apply flag for source image
  void   SetSourceTransformApply(bool);

  /// Get transformation apply flag for source image
  bool GetSourceTransformApply();

  /// Set transformation invert flag for source image
  void   SetSourceTransformInvert(bool);

  /// Get transformation invert flag for source image
  bool GetSourceTransformInvert();

  /// Get a pointer to target image
  irtkImage *GetTarget();

  /// Get a pointer to source image
  irtkImage *GetSource();

  /// Get a pointer to the lookup table of the target image
  irtkLookupTable *GetTargetLookupTable();

  /// Get a pointer to the lookup table of the source image
  irtkLookupTable *GetSourceLookupTable();

  /// Get a pointer to the lookup table of the subtraction of target and source
  irtkLookupTable *GetSubtractionLookupTable();

  /// Get a pointer to the lookup table of the deformation
  irtkLookupTable *GetDeformationLookupTable();

  /// Get transformation
  irtkTransformation *GetTransformation();

  /// Get local transformation
  irtkMultiLevelFreeFormTransformation *GetMFFD();

  /// Reset the display origin to origin of target image
  void Reset();

  /// Get a pointer to segment table
  irtkSegmentTable *GetSegmentTable();

  /// Set interpolation mode for segmentation image
  void   SetSegmentationInterpolationMode(irtkInterpolationMode);

  /// Get interpolation mode for segmentation image
  irtkInterpolationMode GetSegmentationInterpolationMode();

  /// Get a pointer to segmentation image
  irtkGreyImage *GetSegmentation();

  /// Get a pointer to the lookup table of the segmentation image
  irtkLookupTable *GetSegmentationLookupTable();

  /// Get a pointer to the lookup table of the segmentation image
  irtkVoxelContour *GetVoxelContour();

  /// Turns on the segmentation drawing
  void SegmentationLabelsOn();

  /// Turns off the segmentation drawing
  void SegmentationLabelsOff();

  /// Gets the current display value
  int GetDisplaySegmentationLabels();

  /// Turns on the segmentation drawing
  void SegmentationContoursOn();

  /// Turns off the segmentation drawing
  void SegmentationContoursOff();

  /// Gets the current display value
  int GetDisplaySegmentationContours();

  /// Reset the display origin to origin of segmentation image
  void ResetSegmentation();

  /// Set segmentation mode
  void SegmentationMode(int mode);

  /// Set paint brush width
  void SetPaintBrushWidth(int width);

  /// Set minimum region growing threshold
  void SetRegionGrowingThresholdMinimum(int threshold);

  /// Set maximum region growing threshold
  void SetRegionGrowingThresholdMaximum(int threshold);

  /// Returns segmentation mode
  int GetSegmentationMode();

  /// Returns PaintBrushWidth
  int GetPaintBrushWidth();

  /// Get minimum region growing threshold
  int GetRegionGrowingThresholdMinimum();

  /// Get maximum region growing threshold
  int GetRegionGrowingThresholdMaximum();

  /// Clipping of a drawable to the viewport
  void Clip();

  /// Add target landmark
  void AddTargetLandmark(irtkPoint &, char *);

  /// Add source landmark
  void AddSourceLandmark(irtkPoint &, char *);

  /// Delete target landmark
  void DeleteTargetLandmark(int);

  /// Delete source landmark
  void DeleteSourceLandmark(int);

  /// Insert target landmark
  void InsertTargetLandmark(irtkPoint &, int, char *);

  /// Insert source landmark
  void InsertSourceLandmark(irtkPoint &, int, char *);

  /// Get target landmark
  void GetTargetLandmark(irtkPoint &, int, char *);

  /// Get source landmark
  void GetSourceLandmark(irtkPoint &, int, char *);

  /// Put target landmark
  void PutTargetLandmark(irtkPoint, int, char *);

  /// Put source landmark
  void PutSourceLandmark(irtkPoint, int, char *);

  /// Label target landmark
  void LabelTargetLandmark(int, char *);

  /// Label source landmark
  void LabelSourceLandmark(int, char *);

  /// Return number of target landmarks
  int GetNumberOfTargetLandmarks();

  /// Return number of source landmarks
  int GetNumberOfSourceLandmarks();

  /// Approximate transformation using landmark correspondences
  double FitLandmarks();

  /// Callback method for function keys
  void cb_special(int key, int x, int y,
                  int target_delta, int source_delta);

  /// Callback method for function key information
  void cb_special_info();

  /// Callback method for keyboard events
  void cb_keyboard(unsigned char key);

  /// Callback method for keyboard event information
  void cb_keyboard_info();

};

inline void irtkRView::SourceUpdateOn()
{
  _sourceUpdate = true;
}

inline void irtkRView::SegmentationUpdateOn()
{
  _segmentationUpdate = true;
}


inline int irtkRView::GetWidth()
{
  return _screenX;
}

inline int irtkRView::GetHeight()
{
  return _screenY;
}

inline double irtkRView::GetResolution()
{
  return _resolution;
}

inline void irtkRView::SetResolution(double resolution)
{
  _resolution = resolution;
  this->Initialize();
}

inline void irtkRView::SetViewMode(irtkRViewMode value)
{
  _viewMode = value;
}

inline irtkRViewMode irtkRView::GetViewMode()
{
  return _viewMode;
}


inline irtkConfigViewerMode irtkRView::GetConfigMode()
{
  return _configMode;
}

inline void irtkRView::SetConfigMode(irtkConfigViewerMode configMode)
{
  _configMode = configMode;
}

inline void irtkRView::SetViewMix(double value)
{
  _viewMix = value;
}

inline double irtkRView::GetViewMix()
{
  return _viewMix;
}

inline void irtkRView::SetLineThickness(double value)
{
  _LineThickness = value;
}

inline double irtkRView::GetLineThickness()
{
  return _LineThickness;
}

inline void irtkRView::DisplayTargetContoursOn()
{
  _DisplayTargetContour = true;
}

inline void irtkRView::DisplayTargetContoursOff()
{
  _DisplayTargetContour = false;
}

inline int irtkRView::GetDisplayTargetContours()
{
  return _DisplayTargetContour;
}

inline void irtkRView::DisplaySourceContoursOn()
{
  _DisplaySourceContour = true;
}

inline void irtkRView::DisplaySourceContoursOff()
{
  _DisplaySourceContour = false;
}

inline int irtkRView::GetDisplaySourceContours()
{
  return _DisplaySourceContour;
}

inline void irtkRView::SnapToGridOn()
{
  _SnapToGrid = true;

  // Round origin to nearest voxel
  _targetImage->WorldToImage(_origin_x, _origin_y, _origin_z);
  _origin_x = round(_origin_x);
  _origin_y = round(_origin_y);
  _origin_z = round(_origin_z);
  _targetImage->ImageToWorld(_origin_x, _origin_y, _origin_z);

  // Update viewer
  for (int k = 0; k < _NoOfViewers; k++) {
    _targetImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _sourceImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _segmentationImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _selectionImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
  }

  // Update everything else
  _targetUpdate = true;
  _sourceUpdate = true;
  _segmentationUpdate = true;
  _selectionUpdate = true;

}

inline void irtkRView::SnapToGridOff()
{
  _SnapToGrid = false;
}

inline int irtkRView::GetSnapToGrid()
{
  return _SnapToGrid;
}

inline void irtkRView::DisplayCursorOn()
{
  _DisplayCursor = true;
}

inline void irtkRView::DisplayCursorOff()
{
  _DisplayCursor = false;
}

inline irtkCursorMode irtkRView::GetCursorMode()
{
  return _CursorMode;
}

inline void irtkRView::SetCursorMode(irtkCursorMode mode)
{
  _CursorMode = mode;
}

inline int irtkRView::GetDisplayCursor()
{
  return _DisplayCursor;
}

inline void irtkRView::DisplayAxisLabelsOn()
{
  _DisplayAxisLabels = true;
}

inline void irtkRView::DisplayAxisLabelsOff()
{
  _DisplayAxisLabels = false;
}

inline irtkDisplayMode irtkRView::GetDisplayMode()
{
  return _DisplayMode;
}

inline void irtkRView::SetDisplayMode(irtkDisplayMode mode)
{
  _DisplayMode = mode;
}

inline void irtkRView::DisplayDeformationGridOn()
{
  _DisplayDeformationGrid = true;
}

inline void irtkRView::DisplayDeformationGridOff()
{
  _DisplayDeformationGrid = false;
}

inline int irtkRView::GetDisplayDeformationGrid()
{
  return _DisplayDeformationGrid;
}

inline int irtkRView::GetDisplayDeformationGridResolution()
{
  return _DisplayDeformationGridResolution;
}

inline void irtkRView::SetDisplayDeformationGridResolution(int res)
{
  _DisplayDeformationGridResolution = res;
}

inline void irtkRView::DisplayDeformationPointsOn()
{
  _DisplayDeformationPoints = true;
}

inline void irtkRView::DisplayDeformationPointsOff()
{
  _DisplayDeformationPoints = false;
}

inline int irtkRView::GetDisplayDeformationPoints()
{
  return _DisplayDeformationPoints;
}

inline void irtkRView::DisplayDeformationArrowsOn()
{
  _DisplayDeformationArrows = true;
}

inline void irtkRView::DisplayDeformationArrowsOff()
{
  _DisplayDeformationArrows = false;
}

inline int irtkRView::GetDisplayDeformationArrows()
{
  return _DisplayDeformationArrows;
}

inline void irtkRView::DisplayLandmarksOn()
{
  _DisplayLandmarks = true;
}

inline void irtkRView::DisplayLandmarksOff()
{
  _DisplayLandmarks = false;
}

inline int irtkRView::GetDisplayLandmarks()
{
  return _DisplayLandmarks;
}

inline int irtkRView::GetDisplaySegmentationLabels()
{
  return _DisplaySegmentationLabels;
}

inline int irtkRView::GetDisplaySegmentationContours()
{
  return _DisplaySegmentationContours;
}

inline int irtkRView::GetSegmentationMode()
{
  return _SegmentationMode;
}

inline int irtkRView::GetPaintBrushWidth()
{
  return _PaintBrushWidth;
}

inline int irtkRView::GetRegionGrowingThresholdMinimum()
{
  return _RegionGrowingThresholdMin;
}

inline int irtkRView::GetRegionGrowingThresholdMaximum()
{
  return _RegionGrowingThresholdMax;
}

inline void irtkRView::DisplayROIOn()
{
  _DisplayROI = true;
}

inline void irtkRView::DisplayROIOff()
{
  _DisplayROI = false;
}

inline int irtkRView::GetDisplayROI()
{
  return _DisplayROI;
}

inline void irtkRView::TrackTAGOn()
{
  _TrackTAG = true;
}

inline void irtkRView::TrackTAGOff()
{
  _TrackTAG = false;
}

inline int irtkRView::GetTrackTAG()
{
  return _TrackTAG;
}

inline void irtkRView::ViewTAGOn()
{
  _ViewTAG = true;
}

inline void irtkRView::ViewTAGOff()
{
  _ViewTAG = false;
}

inline int irtkRView::GetViewTAG()
{
  return _ViewTAG;
}

#ifdef HAS_VTK

inline vtkPointSet *irtkRView::GetObject(int i)
{
  if ((i < 0) || (i > _NoOfObjects-1)) {
    cerr << "irtkRView::GetObject: Invalid object: " << i << endl;
    return NULL;
  }
  return _Object[i];
}

inline void irtkRView::ObjectMovieOn()
{
    _ObjectMovie = true;
}

inline void irtkRView::ObjectMovieOff()
{
    _ObjectMovie = false;
}

inline int irtkRView::GetObjectMovie()
{
    return _ObjectMovie;
}

inline void irtkRView::DisplayObjectOn()
{
    _DisplayObject = true;
}

inline void irtkRView::DisplayObjectOff()
{
  _DisplayObject = false;
}

inline int irtkRView::GetDisplayObject()
{
  return _DisplayObject;
}

inline void irtkRView::DisplayObjectWarpOn()
{
  _DisplayObjectWarp = true;
}

inline void irtkRView::DisplayObjectWarpOff()
{
  _DisplayObjectWarp = false;
}

inline int irtkRView::GetDisplayObjectWarp()
{
  return _DisplayObjectWarp;
}

inline void irtkRView::DisplayObjectGridOn()
{
  _DisplayObjectGrid = true;
}

inline void irtkRView::DisplayObjectGridOff()
{
  _DisplayObjectGrid = false;
}

inline int irtkRView::GetDisplayObjectGrid()
{
  return _DisplayObjectGrid;
}

#endif

inline irtkRegionGrowingMode irtkRView::GetRegionGrowingMode()
{
  return _regionGrowingMode;
}

inline void irtkRView::SetRegionGrowingMode(irtkRegionGrowingMode mode)
{
  _regionGrowingMode = mode;
}

inline void irtkRView::FlipXOff()
{
  if (_FlipX == true) {
    _FlipX = false;
    _xaxis[0] *= -1;
    _xaxis[1] *= -1;
    _xaxis[2] *= -1;
    this->Initialize();
  }
}

inline void irtkRView::FlipXOn()
{
  if (_FlipX == false) {
    _FlipX = true;
    _xaxis[0] *= -1;
    _xaxis[1] *= -1;
    _xaxis[2] *= -1;
    this->Initialize();
  }
}

inline void irtkRView::FlipYOff()
{
  if (_FlipY == true) {
    _FlipY = false;
    _yaxis[0] *= -1;
    _yaxis[1] *= -1;
    _yaxis[2] *= -1;
    this->Initialize();
  }
}

inline void irtkRView::FlipYOn()
{
  if (_FlipY == false) {
    _FlipY = true;
    _yaxis[0] *= -1;
    _yaxis[1] *= -1;
    _yaxis[2] *= -1;
    this->Initialize();
  }
}

inline void irtkRView::FlipZOff()
{
  if (_FlipZ == true) {
    _FlipZ = false;
    _zaxis[0] *= -1;
    _zaxis[1] *= -1;
    _zaxis[2] *= -1;
    this->Initialize();
  }
}

inline void irtkRView::FlipZOn()
{
  if (_FlipZ == false) {
    _FlipZ = true;
    _zaxis[0] *= -1;
    _zaxis[1] *= -1;
    _zaxis[2] *= -1;
    this->Initialize();
  }
}

inline void irtkRView::SetOrigin(double x, double y, double z)
{
  int i;

  _origin_x = x;
  _origin_y = y;
  _origin_z = z;
  for (i = 0; i < _NoOfViewers; i++) {
    _targetImageOutput[i]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _sourceImageOutput[i]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _segmentationImageOutput[i]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _selectionImageOutput[i]->PutOrigin(_origin_x, _origin_y, _origin_z);
  }
  _targetUpdate = true;
  _sourceUpdate = true;
  _segmentationUpdate = true;
  _selectionUpdate = true;
}

inline void irtkRView::GetOrigin(double &x, double &y, double &z)
{
  x = _origin_x;
  y = _origin_y;
  z = _origin_z;
}

inline void irtkRView::SetROI(double x1, double y1, double z1, double x2, double y2, double z2)
{
  _x1 = x1;
  _y1 = y1;
  _z1 = z1;
  _x2 = x2;
  _y2 = y2;
  _z2 = z2;

}

inline void irtkRView::GetROI(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2)
{
  x1 = _x1;
  y1 = _y1;
  z1 = _z1;
  x2 = _x2;
  y2 = _y2;
  z2 = _z2;
}

inline irtkImage *irtkRView::GetTarget()
{
  return _targetImage;
}

inline irtkImage *irtkRView::GetSource()
{
  return _sourceImage;
}

inline irtkGreyImage *irtkRView::GetSegmentation()
{
  return _segmentationImage;
}

inline irtkVoxelContour *irtkRView::GetVoxelContour()
{
  return &_voxelContour;
}

inline irtkSegmentTable *irtkRView::GetSegmentTable()
{
  return _segmentTable;
}

inline irtkLookupTable *irtkRView::GetTargetLookupTable()
{
  return _targetLookupTable;
}

inline irtkLookupTable *irtkRView::GetSourceLookupTable()
{
  return _sourceLookupTable;
}

inline irtkLookupTable *irtkRView::GetSubtractionLookupTable()
{
  return _subtractionLookupTable;
}

inline irtkTransformation *irtkRView::GetTransformation()
{
  return _sourceTransform;
}

inline irtkMultiLevelFreeFormTransformation *irtkRView::GetMFFD()
{
  irtkMultiLevelFreeFormTransformation *transform = dynamic_cast<irtkMultiLevelFreeFormTransformation *>(_sourceTransform);
  return transform;
}

inline void irtkRView::Clip()
{
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glViewport(0, 0, (GLsizei) _screenX, (GLsizei) _screenY);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0.0, (GLdouble) _screenX, 0.0, (GLdouble) _screenY);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

inline void irtkRView::AddTargetLandmark(irtkPoint &point, char *)
{
  // Add landmark as point, ignoring label for now
      _targetLandmarks.Add(point);
}

inline void irtkRView::AddSourceLandmark(irtkPoint &point, char *)
{
  // Add landmark as point, ignoring label for now	
	  _sourceLandmarks.Add(point);
}

inline void irtkRView::DeleteTargetLandmark(int id)
{
  // Delete landmark from list
  if ((id > 0) && (id <= _targetLandmarks.Size())) {
    irtkPoint p = _targetLandmarks(id-1);
    _targetLandmarks.Del(p);
  }
}

inline void irtkRView::DeleteSourceLandmark(int id)
{
  // Delete landmark from list
  if ((id > 0) && (id <= _sourceLandmarks.Size())) {
    irtkPoint p = _sourceLandmarks(id-1);
    _sourceLandmarks.Del(p);
  }
}

inline void irtkRView::InsertTargetLandmark(irtkPoint &point, int id, char *)
{
  // Insert landmark, ignoring label for now
  if (_targetLandmarks.Size() == 0) {
    _targetLandmarks.Add(point);
  } else if ((id > 0) && (id <= _targetLandmarks.Size())) {
    // Would be nice to call an irtkPointSet::Insert method...
    int i;
    irtkPointSet pset(_targetLandmarks);
    _targetLandmarks.Clear();
    for (i = 1; i < id; i++) {
      _targetLandmarks.Add(pset(i-1));
    }
    _targetLandmarks.Add(point);
    for (i = id; i <= pset.Size(); i++) {
      _targetLandmarks.Add(pset(i-1));
    }
  } else {
    cerr << "irtkRView::InsertTargetLandmark : invalid position " << id << endl;
  }
}

inline void irtkRView::InsertSourceLandmark(irtkPoint &point, int id, char *)
{
  // Insert landmark, ignoring label for now
  if (_sourceLandmarks.Size() == 0) {
    _sourceLandmarks.Add(point);
  } else if ((id > 0) && (id <= _sourceLandmarks.Size())) {
    // Would be nice to call an irtkPointSet::Insert method...
    int i;
    irtkPointSet pset(_sourceLandmarks);
    _sourceLandmarks.Clear();
    for (i = 1; i < id; i++) {
      _sourceLandmarks.Add(pset(i-1));
    }
    _sourceLandmarks.Add(point);
    for (i = id; i <= pset.Size(); i++) {
      _sourceLandmarks.Add(pset(i-1));
    }
  } else {
    cerr << "irtkRView::InsertSourceLandmark : invalid position " << id << endl;
  }
}

inline void irtkRView::LabelTargetLandmark(int, char *)
{
  // So far, labelling of irtkPointSet is not possible
}

inline void irtkRView::LabelSourceLandmark(int, char *)
{
  // So far, labelling of irtkPointSet is not possible
}

inline void irtkRView::GetTargetLandmark(irtkPoint &point, int id, char *)
{
  // Get landmark from list, ignoring label for now
  if ((id > 0) && (id <= _targetLandmarks.Size())) {
    point = _targetLandmarks(id-1);
  }
}

inline void irtkRView::GetSourceLandmark(irtkPoint &point, int id, char *)
{
  // Get landmark from list, ignoring label for now
  if ((id > 0) && (id <= _sourceLandmarks.Size())) {
    point = _sourceLandmarks(id-1);
  }
}

inline void irtkRView::PutTargetLandmark(irtkPoint point, int id, char *)
{
  // Put landmark in list, ignoring label for now
  if ((id > 0) && (id <= _targetLandmarks.Size())) {
    _targetLandmarks(id-1) = point;
  }
}

inline void irtkRView::PutSourceLandmark(irtkPoint point, int id, char *)
{
  // Put landmark in list, ignoring label for now
  if ((id > 0) && (id <= _sourceLandmarks.Size())) {
    _sourceLandmarks(id-1) = point;
  }
}

inline int irtkRView::GetNumberOfTargetLandmarks()
{
  return _targetLandmarks.Size();
}

inline int irtkRView::GetNumberOfSourceLandmarks()
{
  return _sourceLandmarks.Size();
}

inline void irtkRView::SegmentationLabelsOn()
{
  _DisplaySegmentationLabels = true;
}

inline void irtkRView::SegmentationLabelsOff()
{
  _DisplaySegmentationLabels = false;
}

inline void irtkRView::SegmentationContoursOn()
{
  _DisplaySegmentationContours = true;
}

inline void irtkRView::SegmentationContoursOff()
{
  _DisplaySegmentationContours = false;
}

inline double irtkRView::GetTargetMin()
{
  return _targetMin;
}

inline double irtkRView::GetTargetMax()
{
  return _targetMax;
}

inline double irtkRView::GetSourceMin()
{
  return _sourceMin;
}

inline double irtkRView::GetSourceMax()
{
  return _sourceMax;
}

inline double irtkRView::GetSubtractionMin()
{
  return _subtractionMin;
}

inline double irtkRView::GetSubtractionMax()
{
  return _subtractionMax;
}

inline double irtkRView::GetDisplayMinTarget()
{
	return _targetDisplayMin;
}

inline double irtkRView::GetDisplayMaxTarget()
{
	return _targetDisplayMax;
}

inline void irtkRView::SetDisplayMinTarget(double value)
{
	_targetDisplayMin = value;
	_targetLookupTable->SetMinDisplayIntensity(round((value - _targetMin) * 10000.0 / (_targetMax - _targetMin)));
}

inline void irtkRView::SetDisplayMaxTarget(double value)
{
	_targetDisplayMax = value;
	_targetLookupTable->SetMaxDisplayIntensity(round((value - _targetMin) * 10000.0 / (_targetMax - _targetMin)));
}

inline double irtkRView::GetDisplayMinSource()
{
	return _sourceDisplayMin;
}

inline double irtkRView::GetDisplayMaxSource()
{
	return _sourceDisplayMax;
}

inline void irtkRView::SetDisplayMinSource(double value)
{
	_sourceDisplayMin = value;
	_sourceLookupTable->SetMinDisplayIntensity(round((value - _sourceMin) * 10000.0 / (_sourceMax - _sourceMin)));
}

inline void irtkRView::SetDisplayMaxSource(double value)
{
	_sourceDisplayMax = value;
	_sourceLookupTable->SetMaxDisplayIntensity(round((value - _sourceMin) * 10000.0 / (_sourceMax - _sourceMin)));
}

inline double irtkRView::GetDisplayMinSubtraction()
{
	return _subtractionDisplayMin;
}

inline double irtkRView::GetDisplayMaxSubtraction()
{
	return _subtractionDisplayMax;
}

inline void irtkRView::SetDisplayMinSubtraction(double value)
{
	_subtractionDisplayMin = value;
	_subtractionLookupTable->SetMinDisplayIntensity(round(value * 20000.0 / (_subtractionMax - _subtractionMin)));
}

inline void irtkRView::SetDisplayMaxSubtraction(double value)
{
	_subtractionDisplayMax = value;
	_subtractionLookupTable->SetMaxDisplayIntensity(round(value * 20000.0 / (_subtractionMax - _subtractionMin)));
}

#endif
