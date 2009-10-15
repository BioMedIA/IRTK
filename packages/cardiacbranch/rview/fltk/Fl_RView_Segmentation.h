/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

private:

/// Label
char _label[256];

/// ID
int _id;

/// whether the line "_id " is selected
int _selected;

/// Red
unsigned char _red;

/// Green
unsigned char _green;

/// Blue
unsigned char _blue;

/// Transparency
double _trans;

/// Visibility
int _vis;

// main segment panel widgets

/// Menu for paint brush width
static Fl_Menu_Item menu_PaintBrushWidth[];

/// Menu for paint brush width
static Fl_Menu_Item menu_RegionGrowing[];

/// Widget for object browser
Fl_Multi_Browser *segmentObjectBrowser;

///Widget for displaying tha segmentation
Fl_Check_Button  *displayLabels;

/// Widget for the color selected
Fl_Button *editSegmentColor;

/// Widget for label
Fl_Input *editSegmentLabel;

/// Widget for transperancy
Fl_Slider *editSegmentTransperancy;

/// Widget for visibility
Fl_Check_Button *editSegmentVisibility;

/// Window for displaying histograms
Fl_HistogramWindow *_histogramWindow;

/// Button for paintbrush
Fl_Button *drawPaintBrush;

/// Button for contour
Fl_Button *drawContour;

/// Slider for region growing
Fl_Value_Slider2 *regionGrowingMin;

/// Slider for region growing
Fl_Value_Slider2 *regionGrowingMax;

/// Menu for paint brush width
Fl_Choice *paintBrushWidth;

/// Menu for region growing mode
Fl_Choice *regionGrowingMode;

// Callbacks for segmentation
static void cb_loadSegmentation(Fl_Button*, void*);
static void cb_saveSegmentation(Fl_Button*, void*);
static void cb_resetSegmentation(Fl_Button*, void*);
static void cb_browseSegmentObject(Fl_Browser*, void*);
static void cb_addSegment(Fl_Button*, void*);
static void cb_deleteSegment(Fl_Button*, void*);
static void cb_displaySegmentLabels(Fl_Check_Button*, void*);
static void cb_displaySegmentContours(Fl_Check_Button*, void*);
static void cb_deselectAll(Fl_Check_Button* , void* );
static void cb_selectAll(Fl_Check_Button* , void* );
static void cb_loadSegmentTableConfig(Fl_Button*, void*);
static void cb_saveSegmentTableConfig(Fl_Button*, void*);
static void cb_resetSegmentTableConfig(Fl_Button*, void*);
static void cb_editSegmentLabel(Fl_Input*, void*);
static void cb_editSegmentTransperancy(Fl_Slider*, void*);
static void cb_editSegmentVisibility(Fl_Check_Button*, void*);
static void cb_editSegmentPickColor(Fl_Button*, void*);
static void cb_showHistogram(Fl_Button* o, void* v);
static void cb_DrawContour(Fl_Button* o, void* v);
static void cb_DrawPaintBrush(Fl_Button* o, void* v);
static void cb_SetPaintBrushWidth(Fl_Button* o, void* v);
static void cb_SetRegionGrowingMode(Fl_Button* o, void* v);
static void cb_SetRegionGrowingThresholdMinimum(Fl_Value_Slider* o, void* v);
static void cb_SetRegionGrowingThresholdMaximum(Fl_Value_Slider* o, void* v);

//static int fillContour();

public:

/// Show the segmentation control window
void ShowSegmentationControlWindow();

/// Update the segmentation browser
void UpdateSegmentationBrowser();

/// Update the segmentation control window
void UpdateSegmentationControlWindow();

/// Initialize the segmentation control window
void InitializeSegmentationControlWindow();

/// Display segmentation labels
void DisplayLabels();

/// Add label to segmentation
void AddSegment(int label);

