/*=========================================================================
 
 Library   : Image Registration Toolkit (IRTK)
 Module    : $Id$
 Copyright : Imperial College, Department of Computing
 Date      : $Date$
 Version   : $Revision$
 Changes   : $Author$
 
 =========================================================================*/

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkRegistration.h>
#include <irtkPointRegistration.h>

#ifdef __APPLE__
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <irtkRView.h>

#ifdef HAS_VTK
#include <vtkPolyDataReader.h>
#include <vtkStructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#endif

irtkRView::irtkRView(int x, int y)
{
  _screenX = x;
  _screenY = y;

  // Default: Reslice at origin
  _origin_x = 0;
  _origin_z = 0;

  // Default: Resolution is 1 mm
  _resolution = 1;

  // Default: Axis are aligned with Cartesian coordinate system
  _xaxis[0] = 1;
  _xaxis[1] = 0;
  _xaxis[2] = 0;
  _yaxis[0] = 0;
  _yaxis[1] = 1;
  _yaxis[2] = 0;
  _zaxis[0] = 0;
  _zaxis[1] = 0;
  _zaxis[2] = 1;

  // Default: No axis flipping
  _FlipX = false;
  _FlipY = false;
  _FlipZ = false;

  // View mode
  _viewMode = View_A;

  // Default: Viewing mix for shutters is 50%
  _viewMix = 0.5;

  // Default: Interpolation is nearest neighbor
  _targetInterpolator = new irtkNearestNeighborInterpolateImageFunction;
  _sourceInterpolator = new irtkNearestNeighborInterpolateImageFunction;
  _segmentationInterpolator = new irtkNearestNeighborInterpolateImageFunction;
  _selectionInterpolator = new irtkNearestNeighborInterpolateImageFunction;

  // Default time frame
  _targetFrame = 0;
  _sourceFrame = 0;

  // Default: Segmentation
  _DisplaySegmentationLabels = false;
  _DisplaySegmentationContours = false;
  _SegmentationMode = false;
  _PaintBrushWidth = 1;

  // Default: No isolines
  _DisplayTargetContour = false;
  _DisplaySourceContour = false;

  // Default: Line Thickness
  _LineThickness = 2;

  // Default: No ROI
  _DisplayROI = false;

  // Default: No TAG
  _TrackTAG = false;

  // Default: Cursor
  _DisplayCursor = true;

  // Default: Snap to grid
  _SnapToGrid = true;

  // Default: Display mode
  _DisplayMode = Neurological;

  // Default: Deformation properties
  _DeformationProperty = NoneDef;
  _DeformationBlending = 0.5;

  // Default: Axis labels
  _DisplayAxisLabels = true;

  // Default: CrossHair
  _CursorMode = CrossHair;

  // Default: No deformation grid
  _DisplayDeformationGrid = false;
  _DisplayDeformationGridResolution = 0;

  // Default: No deformation points
  _DisplayDeformationPoints = false;

  // Default: No deformation arrows
  _DisplayDeformationArrows = false;

  // Default: Contour is displayed in first viewer
  _contourViewer = -1;
  _contourViewerMode = Viewer_XY;

  // Default: No viewers
  _NoOfViewers = 0;

  // Default: No update needed
  _targetUpdate = false;
  _sourceUpdate = false;

#ifdef HAS_SEGMENTATION_PANEL
  _segmentationUpdate = false;
  _selectionUpdate = false;
#endif

  // Initialize landmark display
  _DisplayLandmarks = false;

#ifdef HAS_VTK
  // Initialize object and display
  _NoOfObjects = 0;
  for (int i = 0; i < MAX_NUMBER_OF_OBJECTS; i++) {
    _Object[i] = NULL;
  }
  _DisplayObject = false;
  _DisplayObjectWarp = false;
  _DisplayObjectGrid = false;
  _ObjectMovie = false;
#endif

  // Allocate memory for source and target image
  _targetImage = new irtkGreyImage;
  _sourceImage = new irtkGreyImage;

  // Allocate memory for segmentation
  _segmentationImage = new irtkGreyImage;

  // Allocate memory for segment Table
  _segmentTable = new irtkSegmentTable();

  // Allocate memory for source and target transformations. Note that in this
  // implementation only the source transformation ever changes. The target
  // transformation should always be an identity transformation.
  _targetTransform = new irtkAffineTransformation;
  _sourceTransform = new irtkAffineTransformation;
  _segmentationTransform = new irtkAffineTransformation;
  _selectionTransform = new irtkAffineTransformation;

  // Flag whether transform shoule be applied
  _sourceTransformApply = true;

  // Flag whether transform should be inverted
  _sourceTransformInvert = false;

  // Initialize min and max values
  _targetMin = 0;
  _targetMax = 1;
  _sourceMin = 0;
  _sourceMax = 1;
  _subtractionMin = 0;
  _subtractionMax = 1;
  _targetDisplayMin = 0;
  _targetDisplayMax = 1;
  _sourceDisplayMin = 0;
  _sourceDisplayMax = 1;
  _subtractionDisplayMin = 0;
  _subtractionDisplayMax = 1;

  // Allocate memory for source and target lookup tables
  _targetLookupTable = new irtkLookupTable;
  _sourceLookupTable = new irtkLookupTable;

  // Allocate memory for subtraction lookup table
  _subtractionLookupTable = new irtkLookupTable;

  // Region growing mode
  _regionGrowingMode = RegionGrowing2D;

  // By default configure rview to start with three orthogonal views
  _configMode = _View_XY_XZ_YZ;
  this->Configure(View_XY_XZ_YZ);
}

irtkRView::~irtkRView()
{
#ifdef HAS_VTK
  for (int i = 0; i < MAX_NUMBER_OF_OBJECTS; i++) {
    if (_Object[i] != NULL) _Object[i]->Delete();
  }
#endif
}

void irtkRView::Update()
{
  int i, j, k, l;
  double blendA, blendB;
  irtkColor *ptr3;
  irtkGreyPixel *ptr1, *ptr2, *ptr4, *ptr5;

  // Check whether target and/or source and/or segmentation need updating
  for (l = 0; l < _NoOfViewers; l++) {
    if ((_targetUpdate == true) && (_targetImage->IsEmpty() != true)) {
      _targetTransformFilter[l]->PutSourcePaddingValue(-1);
      _targetTransformFilter[l]->Run();
    }
    if ((_sourceUpdate == true) && (_sourceImage->IsEmpty() != true)) {
      _sourceTransformFilter[l]->PutSourcePaddingValue(-1);
      _sourceTransformFilter[l]->Run();
    }
    if ((_segmentationUpdate == true)
        && (_segmentationImage->IsEmpty() != true)) {
      _segmentationTransformFilter[l]->Run();
    }
    if ((_selectionUpdate == true)
        && (_voxelContour._raster->IsEmpty() != true)) {
      _selectionTransformFilter[l]->Run();
    }
  }

  // No more updating required
  _targetUpdate = false;
  _sourceUpdate = false;
  _segmentationUpdate = false;
  _selectionUpdate = false;

  // Combine target and source image
  for (k = 0; k < _NoOfViewers; k++) {
    ptr1 = _targetImageOutput[k]->GetPointerToVoxels();
    ptr2 = _sourceImageOutput[k]->GetPointerToVoxels();
    ptr3 = _drawable[k];
    ptr4 = _segmentationImageOutput[k]->GetPointerToVoxels();

    switch (_viewMode) {
      case View_A:
        // Only display the target image
        for (j = 0; j < _viewer[k]->GetHeight(); j++) {
          for (i = 0; i < _viewer[k]->GetWidth(); i++) {
            if (*ptr1 >= 0) {
              *ptr3 = _targetLookupTable->lookupTable[*ptr1];
            } else {
              *ptr3 = irtkColor();
            }
            ptr1++;
            ptr3++;
          }
        }
        break;
      case View_B:
        // Only display the source image
        for (j = 0; j < _viewer[k]->GetHeight(); j++) {
          for (i = 0; i < _viewer[k]->GetWidth(); i++) {
            if (*ptr2 >= 0) {
              *ptr3 = _sourceLookupTable->lookupTable[*ptr2];
            } else {
              *ptr3 = irtkColor();
            }
            ptr2++;
            ptr3++;
          }
        }
        break;
      case View_VShutter:
        // Display target and source images with a vertical shutter
        for (j = 0; j < _viewer[k]->GetHeight(); j++) {
          for (i = 0; i < _viewer[k]->GetWidth(); i++) {
            if (i < _viewMix * _viewer[k]->GetWidth()) {
              if (*ptr1 >= 0) {
                *ptr3 = _targetLookupTable->lookupTable[*ptr1];
              } else {
                *ptr3 = irtkColor();
              }
            } else {
              if (*ptr2 >= _sourceMin) {
                *ptr3 = _sourceLookupTable->lookupTable[*ptr2];
              } else {
                *ptr3 = irtkColor();
                ;
              }
            }
            ptr1++;
            ptr2++;
            ptr3++;
          }
        }
        break;
      case View_HShutter:
        // Display target and source images with a horizontal shutter
        for (j = 0; j < _viewer[k]->GetHeight(); j++) {
          if (j < _viewMix * _viewer[k]->GetHeight()) {
            for (i = 0; i < _viewer[k]->GetWidth(); i++) {
              if (*ptr1 >= 0) {
                *ptr3 = _targetLookupTable->lookupTable[*ptr1];
              } else {
                *ptr3 = irtkColor();
              }
              ptr1++;
              ptr2++;
              ptr3++;
            }
          } else {
            for (i = 0; i < _viewer[k]->GetWidth(); i++) {
              if (*ptr2 >= 0) {
                *ptr3 = _sourceLookupTable->lookupTable[*ptr2];
              } else {
                *ptr3 = irtkColor();
              }
              ptr1++;
              ptr2++;
              ptr3++;
            }
          }
        }
        break;
      case View_Subtraction:
        // Display the subtraction of target and source
        for (j = 0; j < _viewer[k]->GetHeight(); j++) {
          for (i = 0; i < _viewer[k]->GetWidth(); i++) {
            if ((*ptr1 >= 0) && (*ptr2 >= 0)) {
              *ptr3 = _subtractionLookupTable->lookupTable[*ptr1 - *ptr2];
            } else {
              *ptr3 = irtkColor();
            }
            ptr1++;
            ptr2++;
            ptr3++;
          }
        }
        break;
      case View_Checkerboard:
        blendA = _viewMix;
        blendB = 1 - blendA;
        // Display target and source images in a checkerboard fashion
        for (j = 0; j < _viewer[k]->GetHeight(); j++) {
          for (i = 0; i < _viewer[k]->GetWidth(); i++) {
            if ((*ptr1 >= 0) && (*ptr2 >= 0)) {
              ptr3->r = int(blendA * _targetLookupTable->lookupTable[*ptr1].r
                            + blendB * _sourceLookupTable->lookupTable[*ptr2].r);
              ptr3->g = int(blendA * _targetLookupTable->lookupTable[*ptr1].g
                            + blendB * _sourceLookupTable->lookupTable[*ptr2].g);
              ptr3->b = int(blendA * _targetLookupTable->lookupTable[*ptr1].b
                            + blendB * _sourceLookupTable->lookupTable[*ptr2].b);
            } else {
              *ptr3 = irtkColor();
            }
            ptr1++;
            ptr2++;
            ptr3++;
          }
        }
        break;
      case View_AoverB:
        // Display target and source images in a checkerboard fashion
        for (j = 0; j < _viewer[k]->GetHeight(); j++) {
          for (i = 0; i < _viewer[k]->GetWidth(); i++) {
            if ((*ptr1 >= 0) && (*ptr2 >= 0)) {
              ptr3->r = int(_targetLookupTable->lookupTable[*ptr1].a
                            * _targetLookupTable->lookupTable[*ptr1].r + (1
                                - _targetLookupTable->lookupTable[*ptr1].a)
                            * _sourceLookupTable->lookupTable[*ptr2].r);
              ptr3->g = int(_targetLookupTable->lookupTable[*ptr1].a
                            * _targetLookupTable->lookupTable[*ptr1].g + (1
                                - _targetLookupTable->lookupTable[*ptr1].a)
                            * _sourceLookupTable->lookupTable[*ptr2].g);
              ptr3->b = int(_targetLookupTable->lookupTable[*ptr1].a
                            * _targetLookupTable->lookupTable[*ptr1].b + (1
                                - _targetLookupTable->lookupTable[*ptr1].a)
                            * _sourceLookupTable->lookupTable[*ptr2].b);
            } else {
              *ptr3 = irtkColor();
            }
            ptr1++;
            ptr2++;
            ptr3++;
          }
        }
        break;
      case View_BoverA:
        // Display target and source images in a checkerboard fashion
        for (j = 0; j < _viewer[k]->GetHeight(); j++) {
          for (i = 0; i < _viewer[k]->GetWidth(); i++) {
            if ((*ptr1 >= 0) && (*ptr2 >= 0)) {
              ptr3->r = int((1 - _sourceLookupTable->lookupTable[*ptr2].a)
                            * _targetLookupTable->lookupTable[*ptr1].r
                            + _sourceLookupTable->lookupTable[*ptr2].a
                            * _sourceLookupTable->lookupTable[*ptr2].r);
              ptr3->g = int((1 - _sourceLookupTable->lookupTable[*ptr2].a)
                            * _targetLookupTable->lookupTable[*ptr1].g
                            + _sourceLookupTable->lookupTable[*ptr2].a
                            * _sourceLookupTable->lookupTable[*ptr2].g);
              ptr3->b = int((1 - _sourceLookupTable->lookupTable[*ptr2].a)
                            * _targetLookupTable->lookupTable[*ptr1].b
                            + _sourceLookupTable->lookupTable[*ptr2].a
                            * _sourceLookupTable->lookupTable[*ptr2].b);
            } else {
              *ptr3 = irtkColor();
            }
            ptr1++;
            ptr2++;
            ptr3++;
          }
        }
        break;
    }

    if (_DisplaySegmentationLabels == true) {
      ptr3 = _drawable[k];
      // Display segmentation on top of all view modes
      for (j = 0; j < _viewer[k]->GetHeight(); j++) {
        for (i = 0; i < _viewer[k]->GetWidth(); i++) {
          if (*ptr4 >= 0) {
            if (_segmentTable->_entry[*ptr4]._visible == true) {
              blendA = _segmentTable->_entry[*ptr4]._trans;
              blendB = 1 - blendA;
              ptr3->r = int((blendB * ptr3->r) + (blendA
                                                  * _segmentTable->_entry[*ptr4]._color.r));
              ptr3->g = int((blendB * ptr3->g) + (blendA
                                                  * _segmentTable->_entry[*ptr4]._color.g));
              ptr3->b = int((blendB * ptr3->b) + (blendA
                                                  * _segmentTable->_entry[*ptr4]._color.b));
            }
          }
          ptr3++;
          ptr4++;
        }
      }
    }

    if (_voxelContour.Size() > 0) {
      ptr3 = _drawable[k];
      ptr5 = _selectionImageOutput[k]->GetPointerToVoxels();
      // Display segmentation on top of all view modes
      for (j = 0; j < _viewer[k]->GetHeight(); j++) {
        for (i = 0; i < _viewer[k]->GetWidth(); i++) {
          if (*ptr5 > 0) {
            ptr3->r = int((0.5 * ptr3->r) + 0.5 * 255);
            ptr3->g = int((0.5 * ptr3->g) + 0.5 * 255);
            ptr3->b = int((0.5 * ptr3->b));
          }
          ptr3++;
          ptr5++;
        }
      }
    }
  }
}

void irtkRView::Draw()
{
  int k;

  // Clear window
  glClear( GL_COLOR_BUFFER_BIT);

  // Draw images
  for (k = 0; k < _NoOfViewers; k++) {

    // Draw the image
    _viewer[k]->DrawImage(_drawable[k]);

    // Make sure to clip everything to this viewer
    _viewer[k]->Clip();

    // Draw iso-contours in target image if needed
    if (_DisplayTargetContour == true) {
      _viewer[k]->DrawIsolines(_targetImageOutput[k],
                               _targetLookupTable->GetMinDisplayIntensity());
    }
    // Draw iso-contours in source image if needed
    if (_DisplaySourceContour == true) {
      _viewer[k]->DrawIsolines(_sourceImageOutput[k],
                               _sourceLookupTable->GetMinDisplayIntensity());
    }
    // Draw segmentation if needed
    if (_DisplaySegmentationContours == true) {
      _viewer[k]->DrawSegmentationContour(_segmentationImageOutput[k]);
    }
    // Draw tag grid if needed
    if (_ViewTAG == true) {
      // Update grid information based on landmarks
      if(_viewer[k]->UpdateTagGrid(_sourceImageOutput[k], _sourceTransform,_targetLandmarks) == true)
        // If there are 4 landmarks
        _viewer[k]->DrawTagGrid();
    }

    // Update image viewer if necessary
    if ((_DisplayDeformationGrid == true)
        || (_DisplayDeformationPoints == true) || (_DisplayDeformationArrows
            == true)) {
      if (_viewer[k]->Update(_sourceImageOutput[k], _sourceTransform) == true) {

        // Draw deformation grid if needed
        if (_DisplayDeformationGrid == true) {
          _viewer[k]->DrawGrid();
        }
        // Draw deformation points if needed
        if (_DisplayDeformationPoints == true) {
          _viewer[k]->DrawPoints();
        }
        // Draw deformation arrows if needed
        if (_DisplayDeformationArrows == true) {
          _viewer[k]->DrawArrows();
        }
      }
    }

    // Draw landmarks if needed (true: red, false: green)
    if (_DisplayLandmarks == true) {
      _viewer[k]->DrawLandmarks(_targetLandmarks, _targetImageOutput[k], true);
      _viewer[k]->DrawLandmarks(_sourceLandmarks, _targetImageOutput[k], false);
    }

    // Draw ROI if needed
    if (_DisplayROI == true) {
      _viewer[k]->DrawROI(_targetImageOutput[k], _x1, _y1, _z1, _x2, _y2, _z2);
    }

#ifdef HAS_VTK
    // Draw  object if needed
    if (_DisplayObject == true) {
        if(_ObjectMovie == true){
            int _objectFrame = 0;
            if(_targetFrame > _NoOfObjects - 1){
                _objectFrame = _NoOfObjects - 1;
            }else{
                _objectFrame = _targetFrame;
            }
            _viewer[k]->DrawObject(_Object[_objectFrame], _targetImageOutput[k]);
        }else{
          _viewer[k]->DrawObject(_Object, _targetImageOutput[k],
                             _DisplayObjectWarp, _DisplayObjectGrid, _sourceTransform);
        }
    }
#endif

    // Draw cross hairs if needed
    if (_DisplayCursor == true) {
      _viewer[k]->DrawCursor(_CursorMode);
    }

    // Draw axis labels if needed
    if (_DisplayAxisLabels == true) {
      _viewer[k]->DrawInfo(_DisplayMode);
    }

    this->Clip();
  }
}

void irtkRView::SetOrigin(int i, int j)
{
  int k;
  double x1, y1, x2, y2;

  // Convert pixels to normalized coordinates
  _origin_x = i / (double) _screenX;
  _origin_y = (_screenY - j) / (double) _screenY;
  for (k = 0; k < _NoOfViewers; k++) {
    _viewer[k]->GetViewport(x1, y1, x2, y2);
    if ((_origin_x >= x1) && (_origin_x < x2) && (_origin_y >= y1)
        && (_origin_y < y2)) {
      _origin_x = (_origin_x - x1) / (x2 - x1) * _viewer[k]->GetWidth();
      _origin_y = (_origin_y - y1) / (y2 - y1) * _viewer[k]->GetHeight();
      _origin_x = _origin_x;
      _origin_y = _origin_y;
      _origin_z = 0;
      _targetImageOutput[k]->ImageToWorld(_origin_x, _origin_y, _origin_z);
    }
  }

  if (_SnapToGrid == true) {
    // Round origin to nearest voxel
    _targetImage->WorldToImage(_origin_x, _origin_y, _origin_z);
    _origin_x = round(_origin_x);
    _origin_y = round(_origin_y);
    _origin_z = round(_origin_z);
    _targetImage->ImageToWorld(_origin_x, _origin_y, _origin_z);
  }

  for (k = 0; k < _NoOfViewers; k++) {
    _targetImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _sourceImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _segmentationImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _selectionImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
  }

  // Update of target and source is required
  _targetUpdate = true;
  _sourceUpdate = true;
  _segmentationUpdate = true;
  _selectionUpdate = true;

}

void irtkRView::ResetROI()
{
  // Find bounding box
  _x1 = 0;
  _y1 = 0;
  _z1 = 0;
  _targetImage->ImageToWorld(_x1, _y1, _z1);
  _x2 = _targetImage->GetX() - 1;
  _y2 = _targetImage->GetY() - 1;
  _z2 = _targetImage->GetZ() - 1;
  _targetImage->ImageToWorld(_x2, _y2, _z2);
}

void irtkRView::UpdateROI1(int i, int j)
{
  int k;
  double x1, y1, x2, y2, roi1_x, roi1_y, roi1_z, roi2_x, roi2_y, roi2_z;

  // Convert pixels to normalized coordinates
  roi1_x = i / (double) _screenX;
  roi1_y = (_screenY - j) / (double) _screenY;

  // Convert other corner of ROI
  roi2_x = _x2;
  roi2_y = _y2;
  roi2_z = _z2;
  for (k = 0; k < _NoOfViewers; k++) {
    _viewer[k]->GetViewport(x1, y1, x2, y2);
    if ((roi1_x >= x1) && (roi1_x < x2) && (roi1_y >= y1) && (roi1_y < y2)) {
      roi1_x = _x1;
      roi1_y = _y1;
      roi1_z = _z1;
      _targetImageOutput[k]->WorldToImage(roi1_x, roi1_y, roi1_z);
      roi1_x = i / (double) _screenX;
      roi1_y = (_screenY - j) / (double) _screenY;
      roi1_x = (roi1_x - x1) / (x2 - x1) * _viewer[k]->GetWidth();
      roi1_y = (roi1_y - y1) / (y2 - y1) * _viewer[k]->GetHeight();
      _targetImageOutput[k]->ImageToWorld(roi1_x, roi1_y, roi1_z);
      _targetImage->WorldToImage(roi1_x, roi1_y, roi1_z);
      _targetImage->WorldToImage(roi2_x, roi2_y, roi2_z);
      if (round(roi1_x) < 0)
        roi1_x = 0;
      if (round(roi1_x) > round(roi2_x))
        roi1_x = round(roi2_x);
      if (round(roi1_y) < 0)
        roi1_y = 0;
      if (round(roi1_y) > round(roi2_y))
        roi1_y = round(roi2_y);
      if (round(roi1_z) < 0)
        roi1_z = 0;
      if (round(roi1_z) > round(roi2_z))
        roi1_z = round(roi2_z);
      _targetImage->ImageToWorld(roi1_x, roi1_y, roi1_z);
      _x1 = roi1_x;
      _y1 = roi1_y;
      _z1 = roi1_z;
    }
  }
}

void irtkRView::UpdateROI2(int i, int j)
{
  int k;
  double x1, y1, x2, y2, roi1_x, roi1_y, roi1_z, roi2_x, roi2_y, roi2_z;

  // Convert pixels to normalized coordinates
  roi2_x = i / (double) _screenX;
  roi2_y = (_screenY - j) / (double) _screenY;

  // Convert other corner of ROI
  roi1_x = _x1;
  roi1_y = _y1;
  roi1_z = _z1;
  for (k = 0; k < _NoOfViewers; k++) {
    _viewer[k]->GetViewport(x1, y1, x2, y2);
    if ((roi2_x >= x1) && (roi2_x < x2) && (roi2_y >= y1) && (roi2_y < y2)) {
      roi2_x = _x2;
      roi2_y = _y2;
      roi2_z = _z2;
      _targetImageOutput[k]->WorldToImage(roi2_x, roi2_y, roi2_z);
      roi2_x = i / (double) _screenX;
      roi2_y = (_screenY - j) / (double) _screenY;
      roi2_x = (roi2_x - x1) / (x2 - x1) * _viewer[k]->GetWidth();
      roi2_y = (roi2_y - y1) / (y2 - y1) * _viewer[k]->GetHeight();
      _targetImageOutput[k]->ImageToWorld(roi2_x, roi2_y, roi2_z);
      _targetImage->WorldToImage(roi1_x, roi1_y, roi1_z);
      _targetImage->WorldToImage(roi2_x, roi2_y, roi2_z);
      if (round(roi2_x) >= _targetImage->GetX())
        roi2_x = _targetImage->GetX() - 1;
      if (round(roi2_x) < round(roi1_x))
        roi2_x = roi1_x;
      if (round(roi2_y) >= _targetImage->GetY())
        roi2_y = _targetImage->GetY() - 1;
      if (round(roi2_y) < round(roi1_y))
        roi2_y = roi1_y;
      if (round(roi2_z) >= _targetImage->GetZ())
        roi2_z = _targetImage->GetZ() - 1;
      if (round(roi2_z) < round(roi1_z))
        roi2_z = roi1_z;
      _targetImage->ImageToWorld(roi2_x, roi2_y, roi2_z);
      _x2 = roi2_x;
      _y2 = roi2_y;
      _z2 = roi2_z;
    }
  }
}

void irtkRView::AddContour(int i, int j, ContourMode mode)
{
  int k;
  double x1, y1, x2, y2, x, y, z;

  // Convert pixels to normalized coordinates
  x = i / (double) _screenX;
  y = (_screenY - j) / (double) _screenY;

  // If this is the first contour point
  if (_voxelContour.Size() == 0) {
    // Determine in which viewer the contour is to be drawn
    for (k = 0; k < _NoOfViewers; k++) {
      _viewer[k]->GetViewport(x1, y1, x2, y2);
      if ((x >= x1) && (x < x2) && (y >= y1) && (y < y2)) {
        _contourViewer = k;
        _contourViewerMode = _viewer[k]->GetViewerMode();
      }
    }
  } else {
    _viewer[_contourViewer]->GetViewport(x1, y1, x2, y2);
    if ((x < x1) || (x >= x2) || (y < y1) || (y >= y2))
      return;
  }

  // Calculate the coordinates of the contour
  _viewer[_contourViewer]->GetViewport(x1, y1, x2, y2);
  x = (x - x1) / (x2 - x1) * _viewer[_contourViewer]->GetWidth();
  y = (y - y1) / (y2 - y1) * _viewer[_contourViewer]->GetHeight();
  z = 0;

  _targetImageOutput[_contourViewer]->ImageToWorld(x, y, z);

  // Initialise contour if necessary
  if (_voxelContour.Size() == 0) {
    _voxelContour.Initialise(this, _targetImageOutput[_contourViewer]);
  }

  // Add point
  switch (mode) {
    case FirstPoint:
      _voxelContour.AddPointSet(irtkPoint(x, y, z), GetPaintBrushWidth());
      break;
    case NewPoint:
      _voxelContour.AddPoint(irtkPoint(x, y, z), GetPaintBrushWidth());
      break;
    case LastPoint:
      if (_SegmentationMode == 0) {
        _voxelContour.Close(irtkPoint(x, y, z), GetPaintBrushWidth());
      } else {
        _voxelContour.AddPoint(irtkPoint(x, y, z), GetPaintBrushWidth());
      }
      break;
  }

  _selectionUpdate = true;
}

void irtkRView::FillArea(int i, int j)
{
  int k;
  double x1, y1, x2, y2, x, y, z;

  // Convert pixels to normalized coordinates
  x = i / (double) _screenX;
  y = (_screenY - j) / (double) _screenY;

  // If this is the first contour point
  if (_voxelContour.Size() == 0) {
    // Determine in which viewer the contour is to be drawn
    for (k = 0; k < _NoOfViewers; k++) {
      _viewer[k]->GetViewport(x1, y1, x2, y2);
      if ((x >= x1) && (x < x2) && (y >= y1) && (y < y2)) {
        _contourViewer = k;
        _contourViewerMode = _viewer[k]->GetViewerMode();
      }
    }
  } else {
    _viewer[_contourViewer]->GetViewport(x1, y1, x2, y2);
    if ((x < x1) || (x >= x2) || (y < y1) || (y >= y2))
      return;
  }

  // Calculate the coordinates of the contour
  _viewer[_contourViewer]->GetViewport(x1, y1, x2, y2);
  x = (x - x1) / (x2 - x1) * _viewer[_contourViewer]->GetWidth();
  y = (y - y1) / (y2 - y1) * _viewer[_contourViewer]->GetHeight();
  z = 0;
  _targetImageOutput[_contourViewer]->ImageToWorld(x, y, z);

  if (_voxelContour.Size() == 0) {
    _voxelContour.Initialise(this, _targetImageOutput[_contourViewer]);
  }
  _voxelContour.FillArea(irtkPoint(x, y, z));

  _selectionUpdate = true;
}

void irtkRView::RegionGrowContour(int i, int j)
{
  int k;
  double x1, y1, x2, y2, x, y, z;

  // Convert pixels to normalized coordinates
  x = i / (double) _screenX;
  y = (_screenY - j) / (double) _screenY;

  // If this is the first contour point
  if (_voxelContour.Size() == 0) {
    // Determine in which viewer the contour is to be drawn
    for (k = 0; k < _NoOfViewers; k++) {
      _viewer[k]->GetViewport(x1, y1, x2, y2);
      if ((x >= x1) && (x < x2) && (y >= y1) && (y < y2)) {
        _contourViewer = k;
        _contourViewerMode = _viewer[k]->GetViewerMode();
      }
    }
  } else {
    _viewer[_contourViewer]->GetViewport(x1, y1, x2, y2);
    if ((x < x1) || (x >= x2) || (y < y1) || (y >= y2))
      return;
  }

  // Calculate the coordinates of the contour
  _viewer[_contourViewer]->GetViewport(x1, y1, x2, y2);
  x = (x - x1) / (x2 - x1) * _viewer[_contourViewer]->GetWidth();
  y = (y - y1) / (y2 - y1) * _viewer[_contourViewer]->GetHeight();
  z = 0;
  _targetImageOutput[_contourViewer]->ImageToWorld(x, y, z);

  if (_voxelContour.Size() == 0) {
    _voxelContour.Initialise(this, _targetImageOutput[_contourViewer]);
  }
  _voxelContour.RegionGrowing(irtkPoint(x, y, z), _RegionGrowingThresholdMin,
                              _RegionGrowingThresholdMax, _regionGrowingMode);
  _selectionUpdate = true;
}

void irtkRView::UndoContour()
{
  _voxelContour.Undo();
  _selectionUpdate = true;
}

void irtkRView::ClearContour()
{
  _voxelContour.Clear();
  _selectionUpdate = true;
}

void irtkRView::FillContour(int fill, int)
{
  int i, j, k;
  irtkPoint p;

  if (_segmentationImage->IsEmpty() == true) {
    // Create image
    _segmentationImage->Initialize(_targetImage->GetImageAttributes());

    // Fill image with zeros
    irtkGreyPixel *ptr = _segmentationImage->GetPointerToVoxels();
    for (i = 0; i < _segmentationImage->GetNumberOfVoxels(); i++) {
      *ptr = 0;
      ptr++;
    }
  }

  for (k = 0; k < _voxelContour._raster->GetZ(); k++) {
    for (j = 0; j < _voxelContour._raster->GetY(); j++) {
      for (i = 0; i < _voxelContour._raster->GetX(); i++) {
        if (_voxelContour._raster->Get(i, j, k) > 0) {
          p._x = i;
          p._y = j;
          p._z = k;
          _voxelContour._raster->ImageToWorld(p);
          _segmentationImage->WorldToImage(p);
          _segmentationImage->Put(round(p._x), round(p._y), round(p._z), fill);
        }
      }
    }
  }
  _voxelContour.Clear();

  // Update images
  _segmentationUpdate = true;
  _selectionUpdate = true;
}

void irtkRView::Read(char *name)
{
  int ok;
  irtkInterpolationMode interpolation;
  char buffer1[255], *buffer2 = NULL;

  // Open file
  ifstream from(name);
  if (!from) {
    cerr << "irtkRView::Read: Can't open file " << name << "\n";
    exit(1);
  }

  // Read parameters
  while (from.peek() != EOF) {

    ok = false;
    read_line(from, buffer1, buffer2);

    // Config mode (defines NoOfViewers automatically)
    if (strstr(buffer1, "configMode") != NULL) {
      if (strcmp(buffer2, "View_XY") == 0) {
        _configMode = _View_XY;
        ok = true;
      } else if (strcmp(buffer2, "View_XZ") == 0) {
        _configMode = _View_XZ;
        ok = true;
      } else if (strcmp(buffer2, "View_YZ") == 0) {
        _configMode = _View_YZ;
        ok = true;
      } else if (strcmp(buffer2, "View_XY_XZ_v") == 0) {
        _configMode = _View_XY_XZ_v;
        ok = true;
      } else if (strcmp(buffer2, "View_XY_YZ_v") == 0) {
        _configMode = _View_XY_YZ_v;
        ok = true;
      } else if (strcmp(buffer2, "View_XZ_YZ_v") == 0) {
        _configMode = _View_XZ_YZ_v;
        ok = true;
      } else if (strcmp(buffer2, "View_XY_XZ_h") == 0) {
        _configMode = _View_XY_XZ_h;
        ok = true;
      } else if (strcmp(buffer2, "View_XY_YZ_h") == 0) {
        _configMode = _View_XY_YZ_h;
        ok = true;
      } else if (strcmp(buffer2, "View_XZ_YZ_h") == 0) {
        _configMode = _View_XZ_YZ_h;
        ok = true;
      } else if (strcmp(buffer2, "View_XY_XZ_YZ") == 0) {
        _configMode = _View_XY_XZ_YZ;
        ok = true;
      }
    }
    // Width of viewer (in pixels)
    if (strstr(buffer1, "screenX") != NULL) {
      _screenY = atoi(buffer2);
      ok = true;
    }
    // Height of viewer (in pixels)
    if (strstr(buffer1, "screenY") != NULL) {
      _screenY = atoi(buffer2);
      ok = true;
    }
    // Display origin (in mm)
    if (strstr(buffer1, "origin_x") != NULL) {
      _origin_x = atof(buffer2);
      ok = true;
    }
    // Display origin (in mm)
    if (strstr(buffer1, "origin_y") != NULL) {
      _origin_y = atof(buffer2);
      ok = true;
    }
    // Display origin (in mm)
    if (strstr(buffer1, "origin_z") != NULL) {
      _origin_z = atof(buffer2);
      ok = true;
    }
    // Display resolution
    if (strstr(buffer1, "resolution") != NULL) {
      _resolution = atof(buffer2);
      ok = true;
    }

    interpolation = Interpolation_NN;

    // Delete old interpolator
    delete _targetInterpolator;
    // Interpolation mode for target image
    if (strstr(buffer1, "targetInterpolationMode") != NULL) {
      if (strcmp(buffer2, "Interpolation_NN") == 0) {
        interpolation = Interpolation_NN;
        ok = true;
      } else if (strcmp(buffer2, "Interpolation_Linear") == 0) {
        interpolation = Interpolation_Linear;
        ok = true;
      } else if (strcmp(buffer2, "Interpolation_C1Spline") == 0) {
        interpolation = Interpolation_CSpline;
        ok = true;
      } else if (strcmp(buffer2, "Interpolation_BSpline") == 0) {
        interpolation = Interpolation_BSpline;
        ok = true;
      } else if (strcmp(buffer2, "Interpolation_Sinc") == 0) {
        interpolation = Interpolation_Sinc;
        ok = true;
      } else {
        cerr << "irtkRView::Read: Unknown interpolation" << endl;
        exit(1);
      }
    }
    // Create new interpolator
    _targetInterpolator = irtkInterpolateImageFunction::New(interpolation,
                          _targetImage);

    // Delete old interpolator
    delete _sourceInterpolator;
    // Interpolation mode for source image
    if (strstr(buffer1, "sourceInterpolationMode") != NULL) {
      if (strcmp(buffer2, "Interpolation_NN") == 0) {
        interpolation = Interpolation_NN;
        ok = true;
      } else if (strcmp(buffer2, "Interpolation_Linear") == 0) {
        interpolation = Interpolation_Linear;
        ok = true;
      } else if (strcmp(buffer2, "Interpolation_C1Spline") == 0) {
        interpolation = Interpolation_CSpline;
        ok = true;
      } else if (strcmp(buffer2, "Interpolation_BSpline") == 0) {
        interpolation = Interpolation_BSpline;
        ok = true;
      } else if (strcmp(buffer2, "Interpolation_Sinc") == 0) {
        interpolation = Interpolation_Sinc;
        ok = true;
      } else {
        cerr << "irtkRView::Read: Unknown interpolation" << endl;
        exit(1);
      }
    }
    // Create new interpolator
    _sourceInterpolator = irtkInterpolateImageFunction::New(interpolation,
                          _sourceImage);

    // Flag for rview mode
    if (strstr(buffer1, "viewMode") != NULL) {
      if (strcmp(buffer2, "View_A") == 0) {
        _viewMode = View_A;
        ok = true;
      } else if (strcmp(buffer2, "View_B") == 0) {
        _viewMode = View_B;
        ok = true;
      } else if (strcmp(buffer2, "View_Checkerboard") == 0) {
        _viewMode = View_Checkerboard;
        ok = true;
      } else if (strcmp(buffer2, "View_Subtraction") == 0) {
        _viewMode = View_Subtraction;
        ok = true;
      } else if (strcmp(buffer2, "View_HShutter") == 0) {
        _viewMode = View_HShutter;
        ok = true;
      } else if (strcmp(buffer2, "View_VShutter") == 0) {
        _viewMode = View_VShutter;
        ok = true;
      }
    }

    // Display viewing mix in shutter viewing mode
    if (strstr(buffer1, "viewMix") != NULL) {
      _viewMix = atof(buffer2);
      ok = true;
    }
    // Flag for display of isolines from target image
    if (strstr(buffer1, "DisplayTargetContour") != NULL) {
      _DisplayTargetContour = atoi(buffer2);
      ok = true;
    }
    // Flag for display of isolines from source image
    if (strstr(buffer1, "DisplaySourceContour") != NULL) {
      _DisplaySourceContour = atoi(buffer2);
      ok = true;
    }
    // Flag for display of cross hair
    if (strstr(buffer1, "DisplayCursor") != NULL) {
      _DisplayCursor = atoi(buffer2);
      ok = true;
    }
    // Cursor mode
    if (strstr(buffer1, "CursorMode") != NULL) {
      if (strcmp(buffer2, "CrossHair") == 0) {
        _CursorMode = CrossHair;
        ok = true;
      } else if (strcmp(buffer2, "CursorX") == 0) {
        _CursorMode = CursorX;
        ok = true;
      } else if (strcmp(buffer2, "CursorV") == 0) {
        _CursorMode = CursorV;
        ok = true;
      } else if (strcmp(buffer2, "CursorBar") == 0) {
        _CursorMode = CursorBar;
        ok = true;
      } else {
        ok = false;
      }
    }
    // Flag for display of deformation grid
    if (strstr(buffer1, "DisplayDeformationGrid") != NULL) {
      _DisplayDeformationGrid = atoi(buffer2);
      ok = true;
    }
    // Flag for display of deformation points
    if (strstr(buffer1, "DisplayDeformationPoints") != NULL) {
      _DisplayDeformationPoints = atoi(buffer2);
      ok = true;
    }
    // Flag for display of deformation arrows
    if (strstr(buffer1, "DisplayDeformationArrows") != NULL) {
      _DisplayDeformationArrows = atoi(buffer2);
      ok = true;
    }
    // Flag for display of landmarks
    if (strstr(buffer1, "DisplayLandmarks") != NULL) {
      _DisplayLandmarks = atoi(buffer2);
      ok = true;
    }
#ifdef HAS_VTK
    // Flag for display of object
    if (strstr(buffer1, "DisplayObject") != NULL) {
      _DisplayObject = atoi(buffer2);
      ok = true;
    }
    // Flag for warping of object
    if (strstr(buffer1, "DisplayObjectWarp") != NULL) {
      _DisplayObjectWarp = atoi(buffer2);
      ok = true;
    }
    // Flag for display of object grid
    if (strstr(buffer1, "DisplayObjectGrid") != NULL) {
      _DisplayObjectGrid = atoi(buffer2);
      ok = true;
    }
#endif

    // LookupTables - could be replaced by irtkLookupTable stream
    // targetLookupTable
    if (strstr(buffer1, "targetLookupTable_min") != NULL) {
      _targetLookupTable->SetMinDisplayIntensity(atoi(buffer2));
      ok = true;
    }
    if (strstr(buffer1, "targetLookupTable_max") != NULL) {
      _targetLookupTable->SetMaxDisplayIntensity(atoi(buffer2));
      ok = true;
    }
    if (strstr(buffer1, "targetLookupTable_mode") != NULL) {
      if (strcmp(buffer2, "ColorMode_Red") == 0) {
        _targetLookupTable->SetColorModeToRed();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Green") == 0) {
        _targetLookupTable->SetColorModeToGreen();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Blue") == 0) {
        _targetLookupTable->SetColorModeToBlue();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Luminance") == 0) {
        _targetLookupTable->SetColorModeToLuminance();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Rainbow") == 0) {
        _targetLookupTable->SetColorModeToRainbow();
        ok = true;
      }
    }
    // sourceLookupTable
    if (strstr(buffer1, "sourceLookupTable_min") != NULL) {
      _sourceLookupTable->SetMinDisplayIntensity(atoi(buffer2));
      ok = true;
    }
    if (strstr(buffer1, "sourceLookupTable_max") != NULL) {
      _sourceLookupTable->SetMaxDisplayIntensity(atoi(buffer2));
      ok = true;
    }
    if (strstr(buffer1, "sourceLookupTable_mode") != NULL) {
      if (strcmp(buffer2, "ColorMode_Red") == 0) {
        _sourceLookupTable->SetColorModeToRed();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Green") == 0) {
        _sourceLookupTable->SetColorModeToGreen();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Blue") == 0) {
        _sourceLookupTable->SetColorModeToBlue();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Luminance") == 0) {
        _sourceLookupTable->SetColorModeToLuminance();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Rainbow") == 0) {
        _sourceLookupTable->SetColorModeToRainbow();
        ok = true;
      }
    }
    // subtractionLookupTable
    if (strstr(buffer1, "subtractionLookupTable_min") != NULL) {
      _subtractionLookupTable->SetMinDisplayIntensity(atoi(buffer2));
      ok = true;
    }
    if (strstr(buffer1, "subtractionLookupTable_max") != NULL) {
      _subtractionLookupTable->SetMaxDisplayIntensity(atoi(buffer2));
      ok = true;
    }
    if (strstr(buffer1, "subtractionLookupTable_mode") != NULL) {
      if (strcmp(buffer2, "ColorMode_Red") == 0) {
        _subtractionLookupTable->SetColorModeToRed();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Green") == 0) {
        _subtractionLookupTable->SetColorModeToGreen();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Blue") == 0) {
        _subtractionLookupTable->SetColorModeToBlue();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Luminance") == 0) {
        _subtractionLookupTable->SetColorModeToLuminance();
        ok = true;
      } else if (strcmp(buffer2, "ColorMode_Rainbow") == 0) {
        _subtractionLookupTable->SetColorModeToRainbow();
        ok = true;
      }
    }

    // Check if we parse every line
#ifdef DEBUG
    if (ok != true) {
      cerr << "irtkRView::Read() : Ignoring line " << buffer1 << endl;
    }
#endif

  }

  // Configure in the end to take all changed parameters into account
  switch (_configMode) {
    case _View_XY:
      this->Configure(View_XY);
      break;
    case _View_XZ:
      this->Configure(View_XZ);
      break;
    case _View_YZ:
      this->Configure(View_YZ);
      break;
    case _View_XY_XZ_v:
      this->Configure(View_XY_XZ_v);
      break;
    case _View_XY_YZ_v:
      this->Configure(View_XY_YZ_v);
      break;
    case _View_XZ_YZ_v:
      this->Configure(View_XZ_YZ_v);
      break;
    case _View_XY_XZ_h:
      this->Configure(View_XY_XZ_h);
      break;
    case _View_XY_YZ_h:
      this->Configure(View_XY_YZ_h);
      break;
    case _View_XZ_YZ_h:
      this->Configure(View_XZ_YZ_h);
      break;
    case _View_XY_XZ_YZ:
      this->Configure(View_XY_XZ_YZ);
      break;
  }

  // Close file
  from.close();
}

void irtkRView::Write(char *name)
{
  // Open file
  ofstream to(name);
  if (!to) {
    cerr << "irtkRView::Write: Can't open file " << name << "\n";
    exit(1);
  }

  to << "\n#\n# irtkRView configuration\n#\n\n";
  // Write out configuration mode (contains number of viewers)
  switch (_configMode) {
    case _View_XY:
      to << "configMode                        = View_XY\n";
      break;
    case _View_XZ:
      to << "configMode                        = View_XZ\n";
      break;
    case _View_YZ:
      to << "configMode                        = View_YZ\n";
      break;
    case _View_XY_XZ_v:
      to << "configMode                        = View_XY_XZ_v\n";
      break;
    case _View_XY_YZ_v:
      to << "configMode                        = View_XY_YZ_v\n";
      break;
    case _View_XZ_YZ_v:
      to << "configMode                        = View_XZ_YZ_v\n";
      break;
    case _View_XY_XZ_h:
      to << "configMode                        = View_XY_XZ_h\n";
      break;
    case _View_XY_YZ_h:
      to << "configMode                        = View_XY_YZ_h\n";
      break;
    case _View_XZ_YZ_h:
      to << "configMode                        = View_XZ_YZ_h\n";
      break;
    case _View_XY_XZ_YZ:
      to << "configMode                        = View_XY_XZ_YZ\n";
      break;
  }
  // Width of viewer  (in pixels)
  to << "screenX                           = " << _screenX << endl;
  // Height of viewer (in pixels)
  to << "screenY                           = " << _screenY << endl;
  // Display origin (in mm)
  to << "origin_x                          = " << _origin_x << endl;
  // Display origin (in mm)
  to << "origin_y                          = " << _origin_y << endl;
  // Display origin (in mm)
  to << "origin_z                          = " << _origin_z << endl;
  // Display resolution
  to << "resolution                        = " << _resolution << endl;

  // Interpolation mode for target image
  switch (this->GetTargetInterpolationMode()) {
    case Interpolation_NN:
      to << "targetInterpolationMode           = Interpolation_NN\n";
      break;
    case Interpolation_Linear:
      to << "targetInterpolationMode           = Interpolation_Linear\n";
      break;
    case Interpolation_CSpline:
      to << "targetInterpolationMode           = Interpolation_C1Spline\n";
      break;
    case Interpolation_BSpline:
      to << "targetInterpolationMode           = Interpolation_BSpline\n";
      break;
    case Interpolation_Sinc:
      to << "targetInterpolationMode           = Interpolation_Sinc\n";
      break;
    default:
      break;
  }
  // Interpolation mode for source image
  switch (this->GetSourceInterpolationMode()) {
    case Interpolation_NN:
      to << "sourceInterpolationMode           = Interpolation_NN\n";
      break;
    case Interpolation_Linear:
      to << "sourceInterpolationMode           = Interpolation_Linear\n";
      break;
    case Interpolation_CSpline:
      to << "sourceInterpolationMode           = Interpolation_C1Spline\n";
      break;
    case Interpolation_BSpline:
      to << "sourceInterpolationMode           = Interpolation_BSpline\n";
      break;
    case Interpolation_Sinc:
      to << "sourceInterpolationMode           = Interpolation_Sinc\n";
      break;
    default:
      break;
  }

  // Display configuration
  to << "\n#\n# Display configuration\n#\n\n";

  // Flag for rview mode
  switch (_viewMode) {
    case View_A:
      to << "viewMode                          = View_A\n";
      break;
    case View_B:
      to << "viewMode                          = View_B\n";
      break;
    case View_Checkerboard:
      to << "viewMode                          = View_Checkerboard\n";
      break;
    case View_Subtraction:
      to << "viewMode                          = View_Subtraction\n";
      break;
    case View_HShutter:
      to << "viewMode                          = View_HShutter\n";
      break;
    case View_VShutter:
      to << "viewMode                         = View_VShutter\n";
      break;
    default:
      break;
  }
  // Display viewing mix in shutter viewing mode
  to << "viewMix                           = " << _viewMix << endl;
  // Flag for display of isolines from target image
  to << "DisplayTargetContour              = " << _DisplayTargetContour << endl;
  // Flag for display of isolines from source image
  to << "DisplaySourceContour              = " << _DisplaySourceContour << endl;
  // Flag for display of cross hair
  to << "DisplayCursor                     = " << _DisplayCursor << endl;
  // cursor mode
  switch (_CursorMode) {
    case CrossHair:
      to << "CursorMode                        = CrossHair\n";
      break;
    case CursorX:
      to << "CursorMode                        = CursorX\n";
      break;
    case CursorV:
      to << "CursorMode                        = CursorV\n";
      break;
    case CursorBar:
      to << "CursorMode                        = CursorBar\n";
      break;
  }
  // Flag for display of deformation grid
  to << "DisplayDeformationGrid            = " << _DisplayDeformationGrid
  << endl;
  // Flag for display of deformation points
  to << "DisplayDeformationPoints          = " << _DisplayDeformationPoints
  << endl;
  // Flag for display of deformation arrows
  to << "DisplayDeformationArrows          = " << _DisplayDeformationArrows
  << endl;
  // Flag for display of landmarks
  to << "DisplayLandmarks                  = " << _DisplayLandmarks << endl;
#ifdef HAS_VTK
  // Flag for display of object
  to << "DisplayObject                     = " << _DisplayObject << endl;
  // Flag for warping of object
  to << "DisplayObjectWarp                 = " << _DisplayObjectWarp << endl;
  // Flag for display of object grid
  to << "DisplayObjectGrid                 = " << _DisplayObjectGrid << endl;
#endif
  // Lookup tables (could be replaced by irtkLookupTable stream)
  to << "\n#\n# LookupTables\n#\n\n";

  // targetLookupTable
  to << "targetLookupTable_minDisplay      = "
  << _targetLookupTable->GetMinDisplayIntensity() << endl;
  to << "targetLookupTable_maxDisplay      = "
  << _targetLookupTable->GetMaxDisplayIntensity() << endl;
  switch (_targetLookupTable->GetColorMode()) {
    case ColorMode_Red:
      to << "targetLookupTable_mode            = ColorMode_Red\n";
      break;
    case ColorMode_Green:
      to << "targetLookupTable_mode            = ColorMode_Green\n";
      break;
    case ColorMode_Blue:
      to << "targetLookupTable_mode            = ColorMode_Blue\n";
      break;
    case ColorMode_Luminance:
      to << "targetLookupTable_mode            = ColorMode_Luminance\n";
      break;
    case ColorMode_Rainbow:
      to << "targetLookupTable_mode            = ColorMode_Rainbow\n";
      break;
    default:
      break;
  }
  // sourceLookupTable
  to << "sourceLookupTable_minDisplay      = "
  << _sourceLookupTable->GetMinDisplayIntensity() << endl;
  to << "sourceLookupTable_maxDisplay      = "
  << _sourceLookupTable->GetMaxDisplayIntensity() << endl;
  switch (_sourceLookupTable->GetColorMode()) {
    case ColorMode_Red:
      to << "sourceLookupTable_mode            = ColorMode_Red\n";
      break;
    case ColorMode_Green:
      to << "sourceLookupTable_mode            = ColorMode_Green\n";
      break;
    case ColorMode_Blue:
      to << "sourceLookupTable_mode            = ColorMode_Blue\n";
      break;
    case ColorMode_Luminance:
      to << "sourceLookupTable_mode            = ColorMode_Luminance\n";
      break;
    case ColorMode_Rainbow:
      to << "sourceLookupTable_mode            = ColorMode_Rainbow\n";
      break;
    default:
      break;
  }
  // subtractionLookupTable
  to << "subtractionLookupTable_minDisplay = "
  << _subtractionLookupTable->GetMinDisplayIntensity() << endl;
  to << "subtractionLookupTable_maxDisplay = "
  << _subtractionLookupTable->GetMaxDisplayIntensity() << endl;
  switch (_subtractionLookupTable->GetColorMode()) {
    case ColorMode_Red:
      to << "subtractionLookupTable_mode       = ColorMode_Red\n";
      break;
    case ColorMode_Green:
      to << "subtractionLookupTable_mode       = ColorMode_Green\n";
      break;
    case ColorMode_Blue:
      to << "subtractionLookupTable_mode       = ColorMode_Blue\n";
      break;
    case ColorMode_Luminance:
      to << "subtractionLookupTable_mode       = ColorMode_Luminance\n";
      break;
    case ColorMode_Rainbow:
      to << "subtractionLookupTable_mode       = ColorMode_Rainbow\n";
      break;
    default:
      break;
  }

  // Close file
  to.close();
}

void irtkRView::ReadTarget(char *name)
{
  // Read target image
  if (_targetImage != NULL)
    delete _targetImage;
  _targetImage = irtkImage::New(name);

  // Find min and max values and initialize lookup table
  _targetImage->GetMinMaxAsDouble(&_targetMin, &_targetMax);
  _targetLookupTable->Initialize(0, 10000);
  _targetDisplayMin = _targetMin;
  _targetDisplayMax = _targetMax;
  _RegionGrowingThresholdMin = _targetMin;
  _RegionGrowingThresholdMax = _targetMax;

  // Initialize lookup table for subtraction
  _subtractionMin = _targetMin - _sourceMax;
  _subtractionMax = _targetMax - _sourceMin;
  _subtractionDisplayMin = _subtractionMin;
  _subtractionDisplayMax = _subtractionMax;
  _subtractionLookupTable->Initialize(-10000, 10000);

  // Find bounding box
  _x1 = 0;
  _y1 = 0;
  _z1 = 0;
  _targetImage->ImageToWorld(_x1, _y1, _z1);
  _x2 = _targetImage->GetX() - 1;
  _y2 = _targetImage->GetY() - 1;
  _z2 = _targetImage->GetZ() - 1;
  _targetImage->ImageToWorld(_x2, _y2, _z2);

  // Delete contour
  if (_voxelContour.Size() > 0)
    _voxelContour.Clear();

  // Reslice
  this->Reset();
}

void irtkRView::ReadTarget(int argc, char **argv)
{
  irtkImage **nimages;
  int i, n, x, y, z;

  // Determine how many volumes we have
  n = argc;

  // Allocate memory
  nimages = new irtkImage *[n];

  // Read images
  cout << "Reading " << argv[0] << endl;
  nimages[0] = irtkImage::New(argv[0]);

  for (i = 1; i < n; i++) {
    cout << "Reading " << argv[i] << endl;
    nimages[i] = irtkImage::New(argv[i]);
    if (!(nimages[0]->GetImageAttributes() == nimages[i]->GetImageAttributes())) {
      cerr << "Mismatch of image geometry in sequence" << endl;
      nimages[0]->Print();
      nimages[i]->Print();
      exit(1);
    }
  }

  // Delete old image
  if (_targetImage != NULL)
    delete _targetImage;

  // Initialize image attributes
  irtkImageAttributes attr = nimages[0]->GetImageAttributes();
  attr._t = n;
  attr._dt = 1;

  // Allocate new image
  if (dynamic_cast<irtkGenericImage<char> *> (nimages[0]) != NULL) {
    _targetImage = new irtkGenericImage<char> (attr);
  } else if (dynamic_cast<irtkGenericImage<unsigned char> *> (nimages[0])
             != NULL) {
    _targetImage = new irtkGenericImage<unsigned char> (attr);
  } else if (dynamic_cast<irtkGenericImage<short> *> (nimages[0]) != NULL) {
    _targetImage = new irtkGenericImage<short> (attr);
  } else if (dynamic_cast<irtkGenericImage<unsigned short> *> (nimages[0])
             != NULL) {
    _targetImage = new irtkGenericImage<unsigned short> (attr);
  } else if (dynamic_cast<irtkGenericImage<float> *> (nimages[0]) != NULL) {
    _targetImage = new irtkGenericImage<float> (attr);
  } else if (dynamic_cast<irtkGenericImage<double> *> (nimages[0]) != NULL) {
    _targetImage = new irtkGenericImage<double> (attr);
  } else {
    cerr << "irtkRView::ReadTarget: Cannot convert image to desired type"
    << endl;
    exit(1);
  }

  for (i = 0; i < _targetImage->GetT(); i++) {
    for (z = 0; z < _targetImage->GetZ(); z++) {
      for (y = 0; y < _targetImage->GetY(); y++) {
        for (x = 0; x < _targetImage->GetX(); x++) {
          _targetImage->PutAsDouble(x, y, z, i,
                                    nimages[i]->GetAsDouble(x, y, z));
        }
      }
    }
  }
  delete[] nimages;

  // Find min and max values and initialize lookup table
  _targetImage->GetMinMaxAsDouble(&_targetMin, &_targetMax);
  _targetLookupTable->Initialize(0, 10000);
  _targetDisplayMin = _targetMin;
  _targetDisplayMax = _targetMax;
  _RegionGrowingThresholdMin = _targetMin;
  _RegionGrowingThresholdMax = _targetMax;

  // Initialize lookup table for subtraction
  _subtractionMin = _targetMin - _sourceMax;
  _subtractionMax = _targetMax - _sourceMin;
  _subtractionDisplayMin = _subtractionMin;
  _subtractionDisplayMax = _subtractionMax;
  _subtractionLookupTable->Initialize(-10000, 10000);

  // Find bounding box
  _x1 = 0;
  _y1 = 0;
  _z1 = 0;
  _targetImage->ImageToWorld(_x1, _y1, _z1);
  _x2 = _targetImage->GetX() - 1;
  _y2 = _targetImage->GetY() - 1;
  _z2 = _targetImage->GetZ() - 1;
  _targetImage->ImageToWorld(_x2, _y2, _z2);

  // Delete contour
  if (_voxelContour.Size() > 0)
    _voxelContour.Clear();

  // Reslice
  this->Reset();
}

void irtkRView::ReadSource(char *name)
{
  // Read source image
  if (_sourceImage != NULL)
    delete _sourceImage;
  _sourceImage = irtkImage::New(name);

  // Find min and max values and initialize lookup table
  _sourceImage->GetMinMaxAsDouble(&_sourceMin, &_sourceMax);
  _sourceLookupTable->Initialize(0, 10000);
  _sourceDisplayMin = _sourceMin;
  _sourceDisplayMax = _sourceMax;

  // Initialize lookup table for subtraction
  _subtractionMin = _targetMin - _sourceMax;
  _subtractionMax = _targetMax - _sourceMin;
  _subtractionLookupTable->Initialize(_subtractionMin, _subtractionMax);
  _subtractionDisplayMin = _subtractionMin;
  _subtractionDisplayMax = _subtractionMax;
  _subtractionLookupTable->Initialize(-10000, 10000);

  // Update of source is required
  _sourceUpdate = true;

  // Initialize
  this->Initialize();
}

void irtkRView::ReadSource(int argc, char **argv)
{
  irtkImage **nimages;
  int i, n, x, y, z;

  // Determine how many volumes we have
  n = argc;

  // Allocate memory
  nimages = new irtkImage *[n];

  cout << "Reading " << argv[0] << endl;
  nimages[0] = irtkImage::New(argv[0]);

  for (i = 1; i < n; i++) {
    cout << "Reading " << argv[i] << endl;
    nimages[i] = irtkImage::New(argv[i]);
    if (!(nimages[0]->GetImageAttributes() == nimages[i]->GetImageAttributes())) {
      cerr << "Mismatch of image geometry in sequence" << endl;
      exit(1);
    }
  }

  // Delete old image
  if (_sourceImage != NULL)
    delete _sourceImage;

  // Initialize image attributes
  irtkImageAttributes attr = nimages[0]->GetImageAttributes();
  attr._t = n;
  attr._dt = 1;

  // Allocate new image
  if (dynamic_cast<irtkGenericImage<char> *> (nimages[0]) != NULL) {
    _sourceImage = new irtkGenericImage<char> (attr);
  } else if (dynamic_cast<irtkGenericImage<unsigned char> *> (nimages[0])
             != NULL) {
    _sourceImage = new irtkGenericImage<unsigned char> (attr);
  } else if (dynamic_cast<irtkGenericImage<short> *> (nimages[0]) != NULL) {
    _sourceImage = new irtkGenericImage<short> (attr);
  } else if (dynamic_cast<irtkGenericImage<unsigned short> *> (nimages[0])
             != NULL) {
    _sourceImage = new irtkGenericImage<unsigned short> (attr);
  } else if (dynamic_cast<irtkGenericImage<float> *> (nimages[0]) != NULL) {
    _sourceImage = new irtkGenericImage<float> (attr);
  } else if (dynamic_cast<irtkGenericImage<double> *> (nimages[0]) != NULL) {
    _sourceImage = new irtkGenericImage<double> (attr);
  } else {
    cerr << "irtkRView::ReadSource: Cannot convert image to desired type"
    << endl;
    exit(1);
  }

  for (i = 0; i < _sourceImage->GetT(); i++) {
    for (z = 0; z < _sourceImage->GetZ(); z++) {
      for (y = 0; y < _sourceImage->GetY(); y++) {
        for (x = 0; x < _sourceImage->GetX(); x++) {
          _sourceImage->PutAsDouble(x, y, z, i,
                                    nimages[i]->GetAsDouble(x, y, z));
        }
      }
    }
  }
  delete[] nimages;

  // Find min and max values and initialize lookup table
  _sourceImage->GetMinMaxAsDouble(&_sourceMin, &_sourceMax);
  _sourceLookupTable->Initialize(0, 10000);
  _sourceDisplayMin = _sourceMin;
  _sourceDisplayMax = _sourceMax;

  // Initialize lookup table for subtraction
  _subtractionMin = _targetMin - _sourceMax;
  _subtractionMax = _targetMax - _sourceMin;
  _subtractionLookupTable->Initialize(_subtractionMin, _subtractionMax);
  _subtractionDisplayMin = _subtractionMin;
  _subtractionDisplayMax = _subtractionMax;
  _subtractionLookupTable->Initialize(-10000, 10000);

  // Update of source is required
  _sourceUpdate = true;

  // Initialize
  this->Initialize();
}

void irtkRView::ReadSegmentation(char *name)
{
  // Read target image
  _segmentationImage->Read(name);

  // Find bounding box
  _x1 = 0;
  _y1 = 0;
  _z1 = 0;
  _segmentationImage->ImageToWorld(_x1, _y1, _z1);
  _x2 = _segmentationImage->GetX() - 1;
  _y2 = _segmentationImage->GetY() - 1;
  _z2 = _segmentationImage->GetZ() - 1;
  _segmentationImage->ImageToWorld(_x2, _y2, _z2);

  // Update of target is required
  _segmentationUpdate = true;
}

void irtkRView::WriteTarget(char *name)
{
  // Write target image
  _targetImage->Write(name);
}

void irtkRView::WriteSource(char *name)
{
  // Write transformed source image
  if (_sourceImage != NULL) {
    if ((_sourceTransformApply == true) && (_sourceTransform != NULL)) {
      // Allocate the new transformed source and transformation filter
      irtkGreyImage *transformedSource = new irtkGreyImage;
      transformedSource->Initialize(_targetImage->GetImageAttributes());
      irtkImageTransformation *imageTransformation =
        new irtkImageTransformation;
      // Transform source image
      imageTransformation->SetInput(_sourceImage, _sourceTransform);
      imageTransformation->SetOutput(transformedSource);
      imageTransformation->PutInterpolator(_sourceInterpolator);
      imageTransformation->Run();
      // Write transformed image
      transformedSource->Write(name);
      // Be good
      delete transformedSource;
      delete imageTransformation;
    } else {
      _sourceImage->Write(name);
    }
  }
}

void irtkRView::WriteSegmentation(char *name)
{
  // Write target image
  _segmentationImage->Write(name);
}

void irtkRView::ReadTransformation(char *name)
{
  int i;

  // Delete the old transformation
  if (_sourceTransform != NULL)
    delete _sourceTransform;

  // Allocate and read the new transformation
  _sourceTransform = irtkTransformation::New(name);

  // If transformation is rigid convert it to affine
  if (strcmp(_sourceTransform->NameOfClass(), "irtkRigidTransformation") == 0) {
    irtkAffineTransformation *tmpTransform = new irtkAffineTransformation;
    for (i = 0; i < _sourceTransform->NumberOfDOFs(); i++) {
      tmpTransform->Put(i, _sourceTransform->Get(i));
    }
    delete _sourceTransform;
    _sourceTransform = tmpTransform;
  }
  _sourceUpdate = true;

  // Set up the filters
  for (i = 0; i < _NoOfViewers; i++) {
    // Delete the old transformation filter
    delete _sourceTransformFilter[i];

    // Allocate the new transformation filter
    _sourceTransformFilter[i] = irtkImageTransformation::New(_sourceTransform);

    // Set inputs and outputs for the transformation filter
    _sourceTransformFilter[i]->SetInput(_sourceImage);
    _sourceTransformFilter[i]->SetOutput(_sourceImageOutput[i]);
    if (_sourceTransformApply == true) {
      _sourceTransformFilter[i]->SetTransformation(_sourceTransform);
    } else {
      _sourceTransformFilter[i]->SetTransformation(_targetTransform);
    }
    _sourceTransformFilter[i]->PutInterpolator(_sourceInterpolator);
    _sourceTransformFilter[i]->PutSourcePaddingValue(_sourceMin - 1);
    if (_sourceTransformInvert == true) {
      _sourceTransformFilter[i]->InvertOn();
    } else {
      _sourceTransformFilter[i]->InvertOff();
    }
  }
  this->Initialize();
}

void irtkRView::WriteTransformation(char *name)
{
  // Write transformation
  if (_sourceTransform != NULL) {
    _sourceTransform->Write(name);
  }
}

void irtkRView::ReadTargetLandmarks(char *name)
{
  // Read target landmarks
  _targetLandmarks.ReadVTK(name);
}

void irtkRView::ReadSourceLandmarks(char *name)
{
  // Read source landmarks
  _sourceLandmarks.ReadVTK(name);
}

void irtkRView::WriteTargetLandmarks(char *name)
{
  // Write target landmarks
  _targetLandmarks.WriteVTK(name);
}

void irtkRView::WriteSourceLandmarks(char *name)
{
  // Write source landmarks
  _sourceLandmarks.WriteVTK(name);
}

#ifdef HAS_VTK
void irtkRView::ReadObject(const char *name)
{
  if (_NoOfObjects >= MAX_NUMBER_OF_OBJECTS) {
    cerr << "irtkRView::ReadObject(): maximum number of objects reached!\n";
    return;
  }

  // Let vtk do its thing
  vtkPolyDataReader *data_reader = vtkPolyDataReader::New();
  data_reader->SetFileName(name);
  data_reader->Update();
  _Object[_NoOfObjects] = data_reader->GetOutput();
  _Object[_NoOfObjects]->Register(_Object[_NoOfObjects]);

  // Increment objects counter
  _NoOfObjects++;

  // Be good
  data_reader->Delete();
}

void irtkRView::RemoveObject()
{
    int i;

    for (i = 0; i < MAX_NUMBER_OF_OBJECTS; i++) {
        if (_Object[i] != NULL) _Object[i]->Delete();
        _Object[i] = NULL;
    }

    _NoOfObjects = 0;
}

#endif

void irtkRView::Reset()
{
  int iaxis, jaxis, kaxis;
  double xaxis[3], yaxis[3], zaxis[3];

  // Get orientation of target image
  _targetImage->GetOrientation(xaxis, yaxis, zaxis);

  // Compute orientation of axis relative to patient
  _targetImage->Orientation(iaxis, jaxis, kaxis);

  if (_DisplayMode == Native) {
    _xaxis[0] = xaxis[0];
    _xaxis[1] = xaxis[1];
    _xaxis[2] = xaxis[2];
    _yaxis[0] = yaxis[0];
    _yaxis[1] = yaxis[1];
    _yaxis[2] = yaxis[2];
    _zaxis[0] = zaxis[0];
    _zaxis[1] = zaxis[1];
    _zaxis[2] = zaxis[2];
  } else {
    if (_DisplayMode == Neurological) {
      switch (iaxis) {
        case IRTK_L2R:
          _xaxis[0] = -xaxis[0];
          _xaxis[1] = -xaxis[1];
          _xaxis[2] = -xaxis[2];
          break;
        case IRTK_R2L:
          _xaxis[0] = xaxis[0];
          _xaxis[1] = xaxis[1];
          _xaxis[2] = xaxis[2];
          break;
        case IRTK_P2A:
          _yaxis[0] = xaxis[0];
          _yaxis[1] = xaxis[1];
          _yaxis[2] = xaxis[2];
          break;
        case IRTK_A2P:
          _yaxis[0] = -xaxis[0];
          _yaxis[1] = -xaxis[1];
          _yaxis[2] = -xaxis[2];
          break;
        case IRTK_I2S:
          _zaxis[0] = xaxis[0];
          _zaxis[1] = xaxis[1];
          _zaxis[2] = xaxis[2];
          break;
        case IRTK_S2I:
          _zaxis[0] = -xaxis[0];
          _zaxis[1] = -xaxis[1];
          _zaxis[2] = -xaxis[2];
          break;
        default:
          cerr << "irtkRView::ResetTarget: Can't work out x-orientation" << endl;
          break;
      }
      switch (jaxis) {
        case IRTK_L2R:
          _xaxis[0] = -yaxis[0];
          _xaxis[1] = -yaxis[1];
          _xaxis[2] = -yaxis[2];
          break;
        case IRTK_R2L:
          _xaxis[0] = yaxis[0];
          _xaxis[1] = yaxis[1];
          _xaxis[2] = yaxis[2];
          break;
        case IRTK_P2A:
          _yaxis[0] = yaxis[0];
          _yaxis[1] = yaxis[1];
          _yaxis[2] = yaxis[2];
          break;
        case IRTK_A2P:
          _yaxis[0] = -yaxis[0];
          _yaxis[1] = -yaxis[1];
          _yaxis[2] = -yaxis[2];
          break;
        case IRTK_I2S:
          _zaxis[0] = yaxis[0];
          _zaxis[1] = yaxis[1];
          _zaxis[2] = yaxis[2];
          break;
        case IRTK_S2I:
          _zaxis[0] = -yaxis[0];
          _zaxis[1] = -yaxis[1];
          _zaxis[2] = -yaxis[2];
          break;
        default:
          cerr << "irtkRView::ResetTarget: Can't work out y-orientation" << endl;
          break;
      }
      switch (kaxis) {
        case IRTK_L2R:
          _xaxis[0] = -zaxis[0];
          _xaxis[1] = -zaxis[1];
          _xaxis[2] = -zaxis[2];
          break;
        case IRTK_R2L:
          _xaxis[0] = zaxis[0];
          _xaxis[1] = zaxis[1];
          _xaxis[2] = zaxis[2];
          break;
        case IRTK_P2A:
          _yaxis[0] = zaxis[0];
          _yaxis[1] = zaxis[1];
          _yaxis[2] = zaxis[2];
          break;
        case IRTK_A2P:
          _yaxis[0] = -zaxis[0];
          _yaxis[1] = -zaxis[1];
          _yaxis[2] = -zaxis[2];
          break;
        case IRTK_I2S:
          _zaxis[0] = zaxis[0];
          _zaxis[1] = zaxis[1];
          _zaxis[2] = zaxis[2];
          break;
        case IRTK_S2I:
          _zaxis[0] = -zaxis[0];
          _zaxis[1] = -zaxis[1];
          _zaxis[2] = -zaxis[2];
          break;
        default:
          cerr << "irtkRView::ResetTarget: Can't work out z-orientation" << endl;
          break;
      }
    } else {
      switch (iaxis) {
        case IRTK_L2R:
          _xaxis[0] = xaxis[0];
          _xaxis[1] = xaxis[1];
          _xaxis[2] = xaxis[2];
          break;
        case IRTK_R2L:
          _xaxis[0] = -xaxis[0];
          _xaxis[1] = -xaxis[1];
          _xaxis[2] = -xaxis[2];
          break;
        case IRTK_P2A:
          _yaxis[0] = xaxis[0];
          _yaxis[1] = xaxis[1];
          _yaxis[2] = xaxis[2];
          break;
        case IRTK_A2P:
          _yaxis[0] = -xaxis[0];
          _yaxis[1] = -xaxis[1];
          _yaxis[2] = -xaxis[2];
          break;
        case IRTK_I2S:
          _zaxis[0] = xaxis[0];
          _zaxis[1] = xaxis[1];
          _zaxis[2] = xaxis[2];
          break;
        case IRTK_S2I:
          _zaxis[0] = -xaxis[0];
          _zaxis[1] = -xaxis[1];
          _zaxis[2] = -xaxis[2];
          break;
        default:
          cerr << "irtkRView::ResetTarget: Can't work out x-orientation" << endl;
          break;
      }
      switch (jaxis) {
        case IRTK_L2R:
          _xaxis[0] = yaxis[0];
          _xaxis[1] = yaxis[1];
          _xaxis[2] = yaxis[2];
          break;
        case IRTK_R2L:
          _xaxis[0] = -yaxis[0];
          _xaxis[1] = -yaxis[1];
          _xaxis[2] = -yaxis[2];
          break;
        case IRTK_P2A:
          _yaxis[0] = yaxis[0];
          _yaxis[1] = yaxis[1];
          _yaxis[2] = yaxis[2];
          break;
        case IRTK_A2P:
          _yaxis[0] = -yaxis[0];
          _yaxis[1] = -yaxis[1];
          _yaxis[2] = -yaxis[2];
          break;
        case IRTK_I2S:
          _zaxis[0] = yaxis[0];
          _zaxis[1] = yaxis[1];
          _zaxis[2] = yaxis[2];
          break;
        case IRTK_S2I:
          _zaxis[0] = -yaxis[0];
          _zaxis[1] = -yaxis[1];
          _zaxis[2] = -yaxis[2];
          break;
        default:
          cerr << "irtkRView::ResetTarget: Can't work out y-orientation" << endl;
          break;
      }
      switch (kaxis) {
        case IRTK_L2R:
          _xaxis[0] = zaxis[0];
          _xaxis[1] = zaxis[1];
          _xaxis[2] = zaxis[2];
          break;
        case IRTK_R2L:
          _xaxis[0] = -zaxis[0];
          _xaxis[1] = -zaxis[1];
          _xaxis[2] = -zaxis[2];
          break;
        case IRTK_P2A:
          _yaxis[0] = zaxis[0];
          _yaxis[1] = zaxis[1];
          _yaxis[2] = zaxis[2];
          break;
        case IRTK_A2P:
          _yaxis[0] = -zaxis[0];
          _yaxis[1] = -zaxis[1];
          _yaxis[2] = -zaxis[2];
          break;
        case IRTK_I2S:
          _zaxis[0] = zaxis[0];
          _zaxis[1] = zaxis[1];
          _zaxis[2] = zaxis[2];
          break;
        case IRTK_S2I:
          _zaxis[0] = -zaxis[0];
          _zaxis[1] = -zaxis[1];
          _zaxis[2] = -zaxis[2];
          break;
        default:
          cerr << "irtkRView::ResetTarget: Can't work out z-orientation" << endl;
          break;
      }
    }
  }

  // Flip X coordinates if desired
  if (_FlipX == true) {
    _xaxis[0] *= -1;
    _xaxis[1] *= -1;
    _xaxis[2] *= -1;
  }

  // Flip Y coordinates if desired
  if (_FlipY == true) {
    _yaxis[0] *= -1;
    _yaxis[1] *= -1;
    _yaxis[2] *= -1;
  }

  // Flip Z coordinates if desired
  if (_FlipZ == true) {
    _zaxis[0] *= -1;
    _zaxis[1] *= -1;
    _zaxis[2] *= -1;
  }

  // Reslice images at the origin of target image (rounded to nearest voxel)
  _targetImage->GetOrigin(_origin_x, _origin_y, _origin_z);
  _targetImage->WorldToImage(_origin_x, _origin_y, _origin_z);
  _origin_x = round(_origin_x);
  _origin_y = round(_origin_y);
  _origin_z = round(_origin_z);
  _targetImage->ImageToWorld(_origin_x, _origin_y, _origin_z);

  // Initialize viewer
  this->Initialize();

  // Update of target and source is required
  _targetUpdate = true;
  _sourceUpdate = true;
  _segmentationUpdate = true;
  _selectionUpdate = true;
}

void irtkRView::Resize(int w, int h)
{
  int i;

  if ((w != _screenX) || (h != _screenY)) {
    _screenX = w;
    _screenY = h;
    for (i = 0; i < _NoOfViewers; i++) {
      _viewer[i]->SetScreen(_screenX, _screenY);
    }
  }

  _sourceUpdate = true;
  _targetUpdate = true;
  _segmentationUpdate = true;
  _selectionUpdate = true;

  this->Clip();
  this->Initialize();

  // Delete old drawables
  for (i = 0; i < _NoOfViewers; i++) {
    if (_drawable[i] != NULL)
      delete _drawable[i];
  }

  // Allocate new drawables
  for (i = 0; i < _NoOfViewers; i++) {
    int n = _targetImageOutput[i]->GetNumberOfVoxels();
    _drawable[i] = new irtkColor[n];
  }

  this->Update();
}

void irtkRView::Configure(irtkRViewConfig config[])
{
  int i;

  // Delete transformation filter, images and viewers
  for (i = 0; i < _NoOfViewers; i++) {
    delete _targetTransformFilter[i];
    delete _sourceTransformFilter[i];
    delete _segmentationTransformFilter[i];
    delete _selectionTransformFilter[i];
    delete _targetImageOutput[i];
    delete _sourceImageOutput[i];
    delete _segmentationImageOutput[i];
    delete _selectionImageOutput[i];
    delete _viewer[i];
    delete _drawable[i];
  }

  // Delete array for transformation filter, images and viewers
  if (_NoOfViewers > 0) {
    delete[] _targetTransformFilter;
    delete[] _sourceTransformFilter;
    delete[] _segmentationTransformFilter;
    delete[] _selectionTransformFilter;
    delete[] _targetImageOutput;
    delete[] _sourceImageOutput;
    delete[] _segmentationImageOutput;
    delete[] _selectionImageOutput;
    delete[] _viewer;
    delete[] _drawable;
  }

  // Calculate number of viewers
  for (i = 0; config[i].xmin >= 0; i++);
  _NoOfViewers = i;

  // Allocate array for transformation filters
  _targetTransformFilter = new irtkImageTransformation*[_NoOfViewers];
  _sourceTransformFilter = new irtkImageTransformation*[_NoOfViewers];
  _segmentationTransformFilter = new irtkImageTransformation*[_NoOfViewers];
  _selectionTransformFilter = new irtkImageTransformation*[_NoOfViewers];

  // Allocate array for images
  _targetImageOutput = new irtkGreyImage*[_NoOfViewers];
  _sourceImageOutput = new irtkGreyImage*[_NoOfViewers];
  _segmentationImageOutput = new irtkGreyImage*[_NoOfViewers];
  _selectionImageOutput = new irtkGreyImage*[_NoOfViewers];

  // Allocate array for viewers
  _viewer = new irtkViewer*[_NoOfViewers];

  // Allocate array for drawables
  _drawable = new irtkColor*[_NoOfViewers];

  // Configure each viewer
  for (i = 0; i < _NoOfViewers; i++) {

    // Allocate _target viewer
    _viewer[i] = new irtkViewer(this, config[i].mode);

    // Configure _target viewer
    _viewer[i]->SetViewport(config[i].xmin, config[i].ymin, config[i].xmax,
                            config[i].ymax);

    _viewer[i]->SetScreen(_screenX, _screenY);

    _targetImageOutput[i] = new irtkGreyImage;
    _targetTransformFilter[i] = irtkImageTransformation::New(_targetTransform);
    _targetTransformFilter[i]->SetInput(_targetImage);
    _targetTransformFilter[i]->SetOutput(_targetImageOutput[i]);
    _targetTransformFilter[i]->SetTransformation(_targetTransform);
    _targetTransformFilter[i]->PutInterpolator(_targetInterpolator);
    _targetTransformFilter[i]->PutSourcePaddingValue(0);

    _sourceImageOutput[i] = new irtkGreyImage;
    _sourceTransformFilter[i] = irtkImageTransformation::New(_sourceTransform);
    _sourceTransformFilter[i]->SetInput(_sourceImage);
    _sourceTransformFilter[i]->SetOutput(_sourceImageOutput[i]);
    if (_sourceTransformApply == true) {
      _sourceTransformFilter[i]->SetTransformation(_sourceTransform);
    } else {
      _sourceTransformFilter[i]->SetTransformation(_targetTransform);
    }
    _sourceTransformFilter[i]->PutInterpolator(_sourceInterpolator);
    _sourceTransformFilter[i]->PutSourcePaddingValue(_sourceMin - 1);
    if (_sourceTransformInvert == true) {
      _sourceTransformFilter[i]->InvertOn();
    } else {
      _sourceTransformFilter[i]->InvertOff();
    }

    _segmentationImageOutput[i] = new irtkGreyImage;
    _segmentationTransformFilter[i] = irtkImageTransformation::New(
                                        _segmentationTransform);
    _segmentationTransformFilter[i]->SetInput(_segmentationImage);
    _segmentationTransformFilter[i]->SetOutput(_segmentationImageOutput[i]);
    _segmentationTransformFilter[i]->SetTransformation(_segmentationTransform);
    _segmentationTransformFilter[i]->PutInterpolator(_segmentationInterpolator);

    _selectionImageOutput[i] = new irtkGreyImage;
    _selectionTransformFilter[i] = irtkImageTransformation::New(
                                     _selectionTransform);
    _selectionTransformFilter[i]->SetInput(_voxelContour._raster);
    _selectionTransformFilter[i]->SetOutput(_selectionImageOutput[i]);
    _selectionTransformFilter[i]->SetTransformation(_selectionTransform);
    _selectionTransformFilter[i]->PutInterpolator(_selectionInterpolator);
  }
  this->Initialize();

  if (_contourViewer != -1) {
    // Delete contour
    _contourViewer = -1;
    for (i = 0; i < _NoOfViewers; i++) {
      if (_viewer[i]->GetViewerMode() == _contourViewerMode)
        _contourViewer = i;
    }
  }

  for (i = 0; i < _NoOfViewers; i++) {
    // Configure OpenGL
    _drawable[i] = new irtkColor[_targetImageOutput[i]->GetNumberOfVoxels()];
  }

  // Update of target and source is required
  _targetUpdate = true;
  _sourceUpdate = true;
  _segmentationUpdate = true;
  _selectionUpdate = true;
}

void irtkRView::GetInfoText(char *buffer1, char *buffer2, char *buffer3,
                            char *buffer4, char *buffer5)
{
  int i, j, k;
  double u, v, w;

  u = _origin_x;
  v = _origin_y;
  w = _origin_z;
  irtkPoint point(u, v, w);
  _targetImage->WorldToImage(u, v, w);
  i = round(u);
  j = round(v);
  k = round(w);
  if ((i >= 0) && (i < _targetImage->GetX()) && (j >= 0) && (j
      < _targetImage->GetY()) && (k >= 0) && (k < _targetImage->GetZ())
      && (_targetFrame >= 0) && (_targetFrame < _targetImage->GetT())) {
    sprintf(buffer1, "% d % d % d", i, j, k);
    sprintf(buffer2, "% .1f % .1f % .1f", point._x, point._y, point._z);
    sprintf(buffer3, "% .2f", _targetImage->GetAsDouble(i, j, k, _targetFrame));
  } else {
    sprintf(buffer1, " ");
    sprintf(buffer2, " ");
    sprintf(buffer3, " ");
  }
  u = _origin_x;
  v = _origin_y;
  w = _origin_z;
  _sourceTransform->Transform(u, v, w);
  point = irtkPoint(u, v, w);
  _sourceImage->WorldToImage(u, v, w);
  i = round(u);
  j = round(v);
  k = round(w);
  if ((i >= 0) && (i < _sourceImage->GetX()) && (j >= 0) && (j
      < _sourceImage->GetY()) && (k >= 0) && (k < _sourceImage->GetZ())
      && (_sourceFrame >= 0) && (_sourceFrame < _sourceImage->GetT())) {
    sprintf(buffer4, "% .2f", _sourceImage->GetAsDouble(i, j, k, _sourceFrame));
  } else {
    sprintf(buffer4, " ");
  }
  u = _origin_x;
  v = _origin_y;
  w = _origin_z;
  _segmentationTransform->Transform(u, v, w);
  point = irtkPoint(u, v, w);
  _segmentationImage->WorldToImage(u, v, w);
  i = round(u);
  j = round(v);
  k = round(w);
  if ((i >= 0) && (i < _segmentationImage->GetX()) && (j >= 0) && (j
      < _segmentationImage->GetY()) && (k >= 0) && (k
          < _segmentationImage->GetZ())) {
    if ((_segmentationImage->Get(i, j, k) > 0) && (_segmentTable->IsValid(
          _segmentationImage->Get(i, j, k)))) {
      sprintf(buffer5, "%s", _segmentTable->GetLabel(_segmentationImage->Get(i,
              j, k)));
    } else {
      sprintf(buffer5, " ");
    }
  } else {
    sprintf(buffer5, " ");
  }
}

void irtkRView::MouseWheel(int i, int j, int wheel)
{
  int k;
  double u, v, w, x1, y1, x2, y2;

  // Convert pixels to normalized coordinates
  u = i / (double) _screenX;
  v = (_screenY - j) / (double) _screenY;
  w = 0;
  for (k = 0; k < _NoOfViewers; k++) {
    _viewer[k]->GetViewport(x1, y1, x2, y2);
    if ((u >= x1) && (u < x2) && (v >= y1) && (v < y2)) {
      _targetImageOutput[k]->GetOrigin(u, v, w);
      _targetImageOutput[k]->WorldToImage(u, v, w);
      w += wheel;
      _targetImageOutput[k]->ImageToWorld(u, v, w);
    }
  }

  _origin_x = u;
  _origin_y = v;
  _origin_z = w;

  if (_SnapToGrid == true) {
    // Round origin to nearest voxel
    _targetImage->WorldToImage(_origin_x, _origin_y, _origin_z);
    _origin_x = round(_origin_x);
    _origin_y = round(_origin_y);
    _origin_z = round(_origin_z);
    _targetImage->ImageToWorld(_origin_x, _origin_y, _origin_z);
  }

  for (k = 0; k < _NoOfViewers; k++) {
    _targetImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _sourceImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _segmentationImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
    _selectionImageOutput[k]->PutOrigin(_origin_x, _origin_y, _origin_z);
  }

  // Update of target and source is required
  _targetUpdate = true;
  _sourceUpdate = true;
  _segmentationUpdate = true;
  _selectionUpdate = true;

}

void irtkRView::MousePosition(int i, int j)
{
  int k;
  double u, v, w, x1, y1, x2, y2;

  // Convert pixels to normalized coordinates
  u = i / (double) _screenX;
  v = (_screenY - j) / (double) _screenY;
  w = 0;
  for (k = 0; k < _NoOfViewers; k++) {
    _viewer[k]->GetViewport(x1, y1, x2, y2);
    if ((u >= x1) && (u < x2) && (v >= y1) && (v < y2)) {
      u = (u - x1) / (x2 - x1) * _viewer[k]->GetWidth();
      v = (v - y1) / (y2 - y1) * _viewer[k]->GetHeight();
      w = 0;
      _targetImageOutput[k]->ImageToWorld(u, v, w);
      _mouseViewer = k;
    }
  }
  irtkPoint point(u, v, w);
  _targetImage->WorldToImage(u, v, w);
  i = round(u);
  j = round(v);
  k = round(w);
  if ((i >= 0) && (i < _targetImage->GetX()) && (j >= 0) && (j
      < _targetImage->GetY()) && (k >= 0) && (k < _targetImage->GetZ())) {
    _mouseX = i;
    _mouseY = j;
    _mouseZ = k;
    _mouseTargetIntensity = _targetImage->GetAsDouble(i, j, k);
  } else {
    _mouseViewer = -1;
  }
}

void irtkRView::GetTransformationText(list<char *> &text)
{
  int i;
  char *ptr, buffer[256];
  double dx, dy, dz, dt;

  ptr = NULL;
  if (strcmp(_sourceTransform->NameOfClass(), "irtkRigidTransformation") == 0) {
    ptr = strdup("Rigid transformation (6 DOF)");
  } else {
    if (strcmp(_sourceTransform->NameOfClass(), "irtkAffineTransformation") == 0) {
      ptr = strdup("Affine transformation (12 DOF)");
    } else {
      if (strcmp(_sourceTransform->NameOfClass(),
                 "irtkMultiLevelFreeFormTransformation") == 0) {
        ptr = strdup("Affine transformation (12 DOF)");
      } else {
        if (strcmp(_sourceTransform->NameOfClass(),
                   "irtkFluidFreeFormTransformation") == 0) {
          ptr = strdup("Affine transformation (12 DOF)");
        } else {
          if (strcmp(_sourceTransform->NameOfClass(),
                     "irtkMultiLevelFreeFormTransformation4D") == 0) {
            ptr = strdup("Affine transformation (12 DOF)");
          } else {
            ptr = strdup("Unknown transformation type");
          }
        }
      }
    }
  }
  text.push_back(ptr);

  // Convert transformation
  irtkMultiLevelFreeFormTransformation *mffd =
    dynamic_cast<irtkMultiLevelFreeFormTransformation *> (_sourceTransform);

  if (mffd != NULL) {
    for (i = 0; i < mffd->NumberOfLevels(); i++) {
      if (strcmp(mffd->GetLocalTransformation(i)->NameOfClass(),
                 "irtkBSplineFreeFormTransformation4D") == 0) {
        irtkBSplineFreeFormTransformation4D
        *ffd =
          dynamic_cast<irtkBSplineFreeFormTransformation4D *> (mffd->GetLocalTransformation(
                i));
        ffd->GetSpacing(dx, dy, dz, dt);
        sprintf(buffer,
                "4D B-Spline FFD: %d (%.2f mm X %.2f mm X %.2f mm X %.2f ms)",
                ffd->NumberOfDOFs(), dx, dy, dz, dt);
      } else {
        if (strcmp(mffd->GetLocalTransformation(i)->NameOfClass(),
                   "irtkBSplineFreeFormTransformation3D") == 0) {
          irtkBSplineFreeFormTransformation
          *ffd =
            dynamic_cast<irtkBSplineFreeFormTransformation *> (mffd->GetLocalTransformation(
                  i));
          ffd->GetSpacing(dx, dy, dz);
          sprintf(buffer, "3D B-Spline FFD: %d (%.2f mm X %.2f mm X %.2f mm)",
                  ffd->NumberOfDOFs(), dx, dy, dz);
        } else {
          if (strcmp(mffd->GetLocalTransformation(i)->NameOfClass(),
                     "irtkLinearFreeFormTransformation") == 0) {
            irtkLinearFreeFormTransformation
            *ffd =
              dynamic_cast<irtkLinearFreeFormTransformation *> (mffd->GetLocalTransformation(
                    i));
            ffd->GetSpacing(dx, dy, dz);
            sprintf(buffer, "3D Linear FFD: %d (%.2f mm X %.2f mm X %.2f mm)",
                    ffd->NumberOfDOFs(), dx, dy, dz);
          } else {
            if (strcmp(mffd->GetLocalTransformation(i)->NameOfClass(),
                       "irtkEigenFreeFormTransformation") == 0) {
              irtkEigenFreeFormTransformation
              *ffd =
                dynamic_cast<irtkEigenFreeFormTransformation *> (mffd->GetLocalTransformation(
                      i));
              ffd->GetSpacing(dx, dy, dz);
              sprintf(buffer, "3D Eigen FFD: %d (%.2f mm X %.2f mm X %.2f mm)",
                      ffd->NumberOfDOFs(), dx, dy, dz);
            } else {
              cerr
              << "irtkRView::GetTransformationText: Unknown transformation type"
              << endl;
              exit(1);
            }
          }
        }
      }
      ptr = strdup(buffer);
      text.push_back(ptr);
    }
  }
}

void irtkRView::Initialize()
{
  int i;
  irtkImageAttributes attr;

  for (i = 0; i < _NoOfViewers; i++) {
    attr._x = _viewer[i]->GetWidth();
    attr._y = _viewer[i]->GetHeight();
    attr._z = 1;
    attr._xorigin = _origin_x;
    attr._yorigin = _origin_y;
    attr._zorigin = _origin_z;
    switch (_viewer[i]->GetViewerMode()) {
      case Viewer_XY:
        attr._dx = 1.0 / _resolution;
        attr._dy = 1.0 / _resolution;
        attr._dz = _targetImage->GetZSize();
        attr._xaxis[0] = _xaxis[0];
        attr._xaxis[1] = _xaxis[1];
        attr._xaxis[2] = _xaxis[2];
        attr._yaxis[0] = _yaxis[0];
        attr._yaxis[1] = _yaxis[1];
        attr._yaxis[2] = _yaxis[2];
        attr._zaxis[0] = _zaxis[0];
        attr._zaxis[1] = _zaxis[1];
        attr._zaxis[2] = _zaxis[2];
        break;
      case Viewer_XZ:
        attr._dx = 1.0 / _resolution;
        attr._dy = 1.0 / _resolution;
        attr._dz = _targetImage->GetYSize();
        attr._xaxis[0] = _xaxis[0];
        attr._xaxis[1] = _xaxis[1];
        attr._xaxis[2] = _xaxis[2];
        attr._yaxis[0] = _zaxis[0];
        attr._yaxis[1] = _zaxis[1];
        attr._yaxis[2] = _zaxis[2];
        attr._zaxis[0] = _yaxis[0];
        attr._zaxis[1] = _yaxis[1];
        attr._zaxis[2] = _yaxis[2];
        break;
      case Viewer_YZ:
        attr._dx = 1.0 / _resolution;
        attr._dy = 1.0 / _resolution;
        attr._dz = _targetImage->GetZSize();
        attr._xaxis[0] = _yaxis[0];
        attr._xaxis[1] = _yaxis[1];
        attr._xaxis[2] = _yaxis[2];
        attr._yaxis[0] = _zaxis[0];
        attr._yaxis[1] = _zaxis[1];
        attr._yaxis[2] = _zaxis[2];
        attr._zaxis[0] = _xaxis[0];
        attr._zaxis[1] = _xaxis[1];
        attr._zaxis[2] = _xaxis[2];
        break;
      default:
        cerr << "Not a valid viewer mode" << endl;
        exit(1);
        break;
    }
    _targetTransformFilter[i]->SetInput(_targetImage);
    _targetTransformFilter[i]->SetOutput(_targetImageOutput[i]);
    _targetTransformFilter[i]->PutScaleFactorAndOffset(10000.0 / (_targetMax
        - _targetMin), -_targetMin * 10000.0 / (_targetMax - _targetMin));
    _sourceTransformFilter[i]->SetInput(_sourceImage);
    _sourceTransformFilter[i]->SetOutput(_sourceImageOutput[i]);
    _sourceTransformFilter[i]->PutScaleFactorAndOffset(10000.0 / (_sourceMax
        - _sourceMin), -_sourceMin * 10000.0 / (_sourceMax - _sourceMin));
    attr._torigin = _targetFrame;
    _targetImageOutput[i]->Initialize(attr);
    attr._torigin = _sourceFrame;
    _sourceImageOutput[i]->Initialize(attr);
    _segmentationImageOutput[i]->Initialize(attr);
    attr._torigin = 0;
    _selectionImageOutput[i]->Initialize(attr);
  }

  // Update of target and source is required
  _targetUpdate = true;
  _sourceUpdate = true;
  _segmentationUpdate = true;
  _selectionUpdate = true;
}

void irtkRView::SetTargetFrame(int t)
{
  int i;
  double xorigin, yorigin, zorigin, torigin;

  // Update target frame
  _targetFrame = t;

  // Update target origin
  torigin = _targetImage->ImageToTime(t);
  for (i = 0; i < _NoOfViewers; i++) {
    _targetImageOutput[i]->GetOrigin(xorigin, yorigin, zorigin);
    _targetImageOutput[i]->PutOrigin(xorigin, yorigin, zorigin, torigin);
  }

  // Update of target is required
  _targetUpdate = true;
}

int irtkRView::GetTargetFrame()
{
  return _targetFrame;
}

void irtkRView::SetSourceFrame(int t)
{
  int i;
  double xorigin, yorigin, zorigin, torigin;

  if(t >= _sourceImage->GetT()){
      t = _sourceImage->GetT() - 1;
  }

  // Update source frame
  _sourceFrame = t;

  // Update source origin
  torigin = _targetImage->ImageToTime(t);
  for (i = 0; i < _NoOfViewers; i++) {
    _sourceImageOutput[i]->GetOrigin(xorigin, yorigin, zorigin);
    _sourceImageOutput[i]->PutOrigin(xorigin, yorigin, zorigin, torigin);
  }

  // Update of source is required
  _sourceUpdate = true;
}

int irtkRView::GetSourceFrame()
{
  return _sourceFrame;
}

void irtkRView::SetTargetInterpolationMode(irtkInterpolationMode value)
{
  int i;

  delete _targetInterpolator;
  _targetInterpolator = irtkInterpolateImageFunction::New(value, _targetImage);
  for (i = 0; i < _NoOfViewers; i++) {
    _targetTransformFilter[i]->PutInterpolator(_targetInterpolator);
  }
  _targetUpdate = true;
}

irtkInterpolationMode irtkRView::GetTargetInterpolationMode()
{
  if (strstr(_targetInterpolator->NameOfClass(),
             "irtkNearestNeighborInterpolateImageFunction") != NULL) {
    return Interpolation_NN;
  }
  if (strstr(_targetInterpolator->NameOfClass(),
             "irtkLinearInterpolateImageFunction") != NULL) {
    return Interpolation_Linear;
  }
  if (strstr(_targetInterpolator->NameOfClass(),
             "irtkBSplineInterpolateImageFunction") != NULL) {
    return Interpolation_BSpline;
  }
  if (strstr(_targetInterpolator->NameOfClass(),
             "irtkCSplineInterpolateImageFunction") != NULL) {
    return Interpolation_CSpline;
  }
  if (strstr(_targetInterpolator->NameOfClass(),
             "irtkSincInterpolateImageFunction") != NULL) {
    return Interpolation_Sinc;
  }
  return Interpolation_NN;
}

void irtkRView::SetSourceInterpolationMode(irtkInterpolationMode value)
{
  int i;

  delete _sourceInterpolator;
  _sourceInterpolator = irtkInterpolateImageFunction::New(value, _sourceImage);
  for (i = 0; i < _NoOfViewers; i++) {
    _sourceTransformFilter[i]->PutInterpolator(_sourceInterpolator);
  }
  _sourceUpdate = true;
}

irtkInterpolationMode irtkRView::GetSourceInterpolationMode()
{
  if (strstr(_sourceInterpolator->NameOfClass(),
             "irtkNearestNeighborInterpolateImageFunction") != NULL) {
    return Interpolation_NN;
  }
  if (strstr(_sourceInterpolator->NameOfClass(),
             "irtkLinearInterpolateImageFunction") != NULL) {
    return Interpolation_Linear;
  }
  if (strstr(_sourceInterpolator->NameOfClass(),
             "irtkBSplineInterpolateImageFunction") != NULL) {
    return Interpolation_BSpline;
  }
  if (strstr(_sourceInterpolator->NameOfClass(),
             "irtkCSplineInterpolateImageFunction") != NULL) {
    return Interpolation_CSpline;
  }
  if (strstr(_sourceInterpolator->NameOfClass(),
             "irtkSincInterpolateImageFunction") != NULL) {
    return Interpolation_Sinc;
  }
  return Interpolation_NN;
}

void irtkRView::SetSourceTransformInvert(bool value)
{
  int i;

  _sourceTransformInvert = value;
  for (i = 0; i < _NoOfViewers; i++) {
    if (_sourceTransformInvert == true) {
      _sourceTransformFilter[i]->InvertOn();
    } else {
      _sourceTransformFilter[i]->InvertOff();
    }
  }
  _sourceUpdate = true;
}

bool irtkRView::GetSourceTransformInvert()
{
  return _sourceTransformInvert;
}

void irtkRView::SetSourceTransformApply(bool value)
{
  int i;

  _sourceTransformApply = value;
  for (i = 0; i < _NoOfViewers; i++) {
    if (_sourceTransformApply == true) {
      _sourceTransformFilter[i]->SetTransformation(_sourceTransform);
    } else {
      _sourceTransformFilter[i]->SetTransformation(_targetTransform);
    }
  }
  _sourceUpdate = true;
}

bool irtkRView::GetSourceTransformApply()
{
  return _sourceTransformApply;
}

void irtkRView::DrawOffscreen(char *filename)
{
  int i, n, index;
  unsigned char *buffer, *ptr;

  // Allocate memory
  buffer = new unsigned char[_screenX * _screenY * 3];

  // Make sure everything is setup correctly (this may be the first time
  // something is drawn into the window)
  this->Resize(_screenX, _screenY);
  this->Draw();

  // Force framebuffer to flush
  glFlush();

  // Allocate RGB image
  irtkGenericImage<unsigned char> image(this->GetWidth(), this->GetHeight(), 3, 1);

  // Read pixels from framebuffer
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, this->GetWidth(), this->GetHeight(), GL_RGB,
               GL_UNSIGNED_BYTE, buffer);

  // Convert to RGB image
  n = image.GetX() * image.GetY();
  index = 0;
  ptr = image.GetPointerToVoxels();
  for (i = 0; i < n; i++) {
    ptr[i]     = buffer[index++];
    ptr[i+n]   = buffer[index++];
    ptr[i+2*n] = buffer[index++];
  }

  // Write file
  image.Write(filename);

  delete[] buffer;
}

double irtkRView::FitLandmarks()
{
  int i, n;
  double error;
  irtkPointSet target_pts, source_pts;

  // Check whether landmarks numbers agree
  if (this->GetNumberOfTargetLandmarks() != this->GetNumberOfSourceLandmarks()) {
    return 0;
  }

  // Initialize
  n = this->GetNumberOfTargetLandmarks();
  for (i = 0; i < n; i++) {
    target_pts.Add(this->_targetLandmarks(i));
    source_pts.Add(this->_sourceLandmarks(i));
  }

  // Set transformation
  irtkRigidTransformation transformation;

#ifdef HAS_VTK
  // Set input and output for the registration filter
  irtkPointRigidRegistration registration;
  registration.SetInput (&target_pts, &source_pts);
  registration.SetOutput(&transformation);

  // Run registration filter
  registration.Run();
#endif

  // Calculate residual error
  transformation.irtkTransformation::Transform(source_pts);

  error = 0;
  for (i = 0; i < target_pts.Size(); i++) {
    irtkPoint p1 = target_pts(i);
    irtkPoint p2 = source_pts(i);
    error += sqrt(pow(p1._x - p2._x, 2.0) + pow(p1._y - p2._y, 2.0) + pow(p1._z
                  - p2._z, 2.0));
  }
  error /= target_pts.Size();

  // Copy dofs
  for (i = 0; i < 6; i++) {
    _sourceTransform->Put(i, transformation.Get(i));
  }

  return error;
}

void irtkRView::cb_special(int key, int, int, int target_delta,
                           int source_delta)
{
  switch (key) {
    case KEY_F1:
      _targetLookupTable->SetMinDisplayIntensity(
        _targetLookupTable->GetMinDisplayIntensity() + target_delta);
      break;
    case KEY_F2:
      _targetLookupTable->SetMinDisplayIntensity(
        _targetLookupTable->GetMinDisplayIntensity() - target_delta);
      break;
    case KEY_F3:
      _targetLookupTable->SetMaxDisplayIntensity(
        _targetLookupTable->GetMaxDisplayIntensity() + target_delta);
      break;
    case KEY_F4:
      _targetLookupTable->SetMaxDisplayIntensity(
        _targetLookupTable->GetMaxDisplayIntensity() - target_delta);
      break;
    case KEY_F5:
      _sourceLookupTable->SetMinDisplayIntensity(
        _sourceLookupTable->GetMinDisplayIntensity() + source_delta);
      break;
    case KEY_F6:
      _sourceLookupTable->SetMinDisplayIntensity(
        _sourceLookupTable->GetMinDisplayIntensity() - source_delta);
      break;
    case KEY_F7:
      _sourceLookupTable->SetMaxDisplayIntensity(
        _sourceLookupTable->GetMaxDisplayIntensity() + source_delta);
      break;
    case KEY_F8:
      _sourceLookupTable->SetMaxDisplayIntensity(
        _sourceLookupTable->GetMaxDisplayIntensity() - source_delta);
      break;
    case KEY_F9:
      if (this->GetDisplayTargetContours()) {
        this->DisplayTargetContoursOff();
      } else {
        this->DisplayTargetContoursOn();
      }
      break;
    case KEY_F10:
      if (this->GetDisplaySourceContours()) {
        this->DisplaySourceContoursOff();
      } else {
        this->DisplaySourceContoursOn();
      }
      break;
    case KEY_F11:
      break;
    case KEY_F12:
      break;
    default:
      break;
  }
  this->Update();
  this->Draw();
}

void irtkRView::cb_special_info()
{
  cerr << "\tSpecial function keys:\n";
  cerr << "\tF1                               Increase min target intensity\n";
  cerr << "\tF2                               Decrease min target intensity\n";
  cerr << "\tF3                               Increase max target intensity\n";
  cerr << "\tF4                               Decrease max target intensity\n";
  cerr << "\tF5                               Increase min source intensity\n";
  cerr << "\tF6                               Decrease min source intensity\n";
  cerr << "\tF7                               Increase max source intensity\n";
  cerr << "\tF8                               Decrease max source intensity\n";
  cerr << "\tF9                               Display target contours on\n";
  cerr << "\tF10                              Display target contours off\n";
  cerr << "\tF11                              Display source contours on\n";
  cerr << "\tF12                              Display source contours off\n\n";
  return;
}

void irtkRView::cb_keyboard(unsigned char key)
{
  double x, y, z, t;

  switch (key) {
    case 27:
      exit(0);
      break;
    case 'q':
      exit(0);
      break;
    case 'r':
      this->Reset();
      break;
    case 'i':
      this->ResetROI();
      break;
    case 'l':
      this->SetTargetInterpolationMode(Interpolation_Linear);
      this->SetSourceInterpolationMode(Interpolation_Linear);
      break;
    case 'n':
      this->SetTargetInterpolationMode(Interpolation_NN);
      this->SetSourceInterpolationMode(Interpolation_NN);
      break;
    case 'c':
      this->SetTargetInterpolationMode(Interpolation_CSpline);
      this->SetSourceInterpolationMode(Interpolation_CSpline);
      break;
    case 'b':
      this->SetTargetInterpolationMode(Interpolation_BSpline);
      this->SetSourceInterpolationMode(Interpolation_BSpline);
      break;
    case 'S':
      this->SetTargetInterpolationMode(Interpolation_Sinc);
      this->SetSourceInterpolationMode(Interpolation_Sinc);
      break;
    case 't':
      this->SetViewMode(View_A);
      break;
    case 's':
      this->SetViewMode(View_B);
      break;
    case 'm':
      this->SetViewMode(View_Checkerboard);
      break;
    case 'd':
      this->SetViewMode(View_Subtraction);
      break;
    case ' ':
      if (this->GetDisplayCursor()) {
        this->DisplayCursorOff();
      } else {
        this->DisplayCursorOn();
      }
      break;
    case 'h':
      this->SetCursorMode(CrossHair);
      break;
    case 'x':
      this->GetOrigin(x, y, z);
      _targetImage->WorldToImage(x, y, z);
      x--;
      if (x < 0)
        x = 0;
      if (x >= _targetImage->GetX())
        x = _targetImage->GetX() - 1;
      _targetImage->ImageToWorld(x, y, z);
      this->SetOrigin(x, y, z);
      break;
    case 'X':
      this->GetOrigin(x, y, z);
      _targetImage->WorldToImage(x, y, z);
      x++;
      if (x < 0)
        x = 0;
      if (x >= _targetImage->GetX())
        x = _targetImage->GetX() - 1;
      _targetImage->ImageToWorld(x, y, z);
      this->SetOrigin(x, y, z);
      break;
    case 'y':
      this->GetOrigin(x, y, z);
      _targetImage->WorldToImage(x, y, z);
      y--;
      if (y < 0)
        y = 0;
      if (y >= _targetImage->GetX())
        y = _targetImage->GetY() - 1;
      _targetImage->ImageToWorld(x, y, z);
      this->SetOrigin(x, y, z);
      break;
    case 'Y':
      this->GetOrigin(x, y, z);
      _targetImage->WorldToImage(x, y, z);
      y++;
      if (y < 0)
        y = 0;
      if (y >= _targetImage->GetX())
        y = _targetImage->GetY() - 1;
      _targetImage->ImageToWorld(x, y, z);
      this->SetOrigin(x, y, z);
      break;
    case 'z':
      this->GetOrigin(x, y, z);
      _targetImage->WorldToImage(x, y, z);
      z--;
      if (z < 0)
        z = 0;
      if (z >= _targetImage->GetX())
        z = _targetImage->GetZ() - 1;
      _targetImage->ImageToWorld(x, y, z);
      this->SetOrigin(x, y, z);
      break;
    case 'Z':
      this->GetOrigin(x, y, z);
      _targetImage->WorldToImage(x, y, z);
      z++;
      if (z < 0)
        z = 0;
      if (z >= _targetImage->GetX())
        z = _targetImage->GetZ() - 1;
      _targetImage->ImageToWorld(x, y, z);
      this->SetOrigin(x, y, z);
      break;
    case '+':
      t = this->GetTargetFrame();
      t++;
      if (t < 0)
        t = _targetImage->GetT() - 1;
      if (t >= _targetImage->GetT())
        t = 0;
      this->SetTargetFrame(round(t));
      break;
    case '-':
      t = this->GetTargetFrame();
      t--;
      if (t < 0)
        t = _targetImage->GetT() - 1;
      if (t >= _targetImage->GetT())
        t = 0;
      this->SetTargetFrame(round(t));
      break;
    case 'v':
      this->SetCursorMode(CursorV);
      break;
    case 'B':
      this->SetCursorMode(CursorBar);
      break;
    case 'g':
      if (this->GetDisplayDeformationGrid()) {
        this->DisplayDeformationGridOff();
      } else {
        this->DisplayDeformationGridOn();
      }
      break;
    case 'p':
      if (this->GetDisplayDeformationPoints()) {
        this->DisplayDeformationPointsOff();
      } else {
        this->DisplayDeformationPointsOn();
      }
      break;
    case 'a':
      if (this->GetDisplayDeformationArrows()) {
        this->DisplayDeformationArrowsOff();
      } else {
        this->DisplayDeformationArrowsOn();
      }
      break;
    case 'L':
      if (this->GetDisplayLandmarks()) {
        this->DisplayLandmarksOff();
      } else {
        this->DisplayLandmarksOn();
      }
      break;
#ifdef HAS_VTK
    case 'O':
      if (this->GetDisplayObject()) {
        this->DisplayObjectOff();
      } else {
        this->DisplayObjectOn();
      }
      break;
    case 'W':
      if (this->GetDisplayObjectWarp()) {
        this->DisplayObjectWarpOff();
      } else {
        this->DisplayObjectWarpOn();
      }
      break;
    case 'G':
      if (this->GetDisplayObjectGrid()) {
        this->DisplayObjectGridOff();
      } else {
        this->DisplayObjectGridOn();
      }
      break;
#endif
    case '>':
      this->SetResolution(this->GetResolution() * 2.0);
      break;
    case '<':
      this->SetResolution(this->GetResolution() / 2.0);
      break;
    case '.':
      this->SetResolution(this->GetResolution() * sqrt(2.0));
      break;
    case ',':
      this->SetResolution(this->GetResolution() / sqrt(2.0));
      break;
    default:
      break;
  }
  this->Update();
  this->Draw();

  return;
}

void irtkRView::SegmentationMode(int mode)
{
  _SegmentationMode = mode;
}

void irtkRView::SetPaintBrushWidth(int width)
{
  _PaintBrushWidth = width;
}

void irtkRView::SetRegionGrowingThresholdMinimum(int threshold)
{
  _RegionGrowingThresholdMin = threshold;
}

void irtkRView::SetRegionGrowingThresholdMaximum(int threshold)
{
  _RegionGrowingThresholdMax = threshold;
}

void irtkRView::cb_keyboard_info()
{
  cerr << "\tControl keys:\n";
  cerr << "\t'q'                              Exit\n";
  cerr << "\t'r'                              Reset target\n";
  cerr << "\t'R'                              Reset source\n";
  cerr << "\t'l'                              Linear interpolation\n";
  cerr
  << "\t'n'                              Nearest neighbour interpolation\n";
  cerr << "\t'c'                              C1-spline interpolation\n";
  cerr << "\t'b'                              B-spline interpolation\n";
  cerr << "\t'S'                              Sinc interpolation\n";
  cerr << "\t't'                              View target\n";
  cerr << "\t's'                              View source\n";
  cerr << "\t'm'                              Mixed viewport (checkerboard)\n";
  cerr << "\t'd'                              View difference (subtraction)\n";
  cerr << "\t' '                              Cursor on/off\n";
  cerr << "\t'h'                              Display cursor as cross hair\n";
  cerr << "\t'x'                              Display cursor as X\n";
  cerr << "\t'v'                              Display cursor as V\n";
  cerr << "\t'B'                              Display cursor as bar\n";
  cerr << "\t'g'                              Deformation grid     on/off\n";
  cerr << "\t'p'                              Deformation points   on/off\n";
#ifndef HAS_VTK
  cerr << "\t'@'                              Deformation labels   on/off\n";
#endif
  cerr << "\t'a'                              Deformation arrows   on/off\n";
  cerr << "\t'='                              Relative deformation on/off\n";
  cerr << "\t'+'                              Increase deformation level\n";
  cerr << "\t'-'                              Decrease deformation level\n";
  cerr << "\t'L'                              Landmarks on/off\n";
#ifdef HAS_VTK
  cerr << "\t'O'                              Object display on/off\n";
  cerr << "\t'W'                              Object vectors warp on/off\n";
  cerr << "\t'G'                              Object grid on/off\n";
#endif
  cerr
  << "\t'>'                              Increase resolution by factor 2\n";
  cerr
  << "\t'<'                              Decrease resolution by factor 1/2\n";
  cerr
  << "\t'.'                              Increase resolution by factor sqrt(2)\n";
  cerr
  << "\t','                              Decrease resolution by factor 1/sqrt(2)\n\n";

  return;
}
