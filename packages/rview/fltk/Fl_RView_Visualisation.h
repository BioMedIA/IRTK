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

//
// Members for image display
//

/// Interpolation menu for target image
static Fl_Menu_Item menu_targetInterpolationMode[];

/// Interpolation menu for source image
static Fl_Menu_Item menu_sourceInterpolationMode[];

/// Interpolation menu for target image
static Fl_Menu_Item menu_targetColorMode[];

/// Interpolation menu for source image
static Fl_Menu_Item menu_sourceColorMode[];

/// Widget for minimum intensity of target image
Fl_Value_Slider2 *targetMin;

/// Widget for maximum intensity of target image
Fl_Value_Slider2 *targetMax;

/// Widget for target image frame
Fl_Value_Slider *targetFrame;

/// Widget for minimum intensity of source image
Fl_Value_Slider2 *sourceMin;

/// Widget for maximum intensity of source image
Fl_Value_Slider2 *sourceMax;

/// Widget for source image frame
Fl_Value_Slider *sourceFrame;

/// Widget for source image frame
Fl_Value_Slider *lineThickness;

/// Widget for movie speed
Fl_Value_Slider *speed;

/// Widget for target interpolation
Fl_Choice *targetInterpolationMode;

/// Widget for source interpolation
Fl_Choice *sourceInterpolationMode;

/// Widget for target color
Fl_Choice *targetColorMode;

/// Widget for source color
Fl_Choice *sourceColorMode;

/// Widget for display control of target isolines
Fl_Check_Button *TargetIsolines;

/// Widget for display control of source isolines
Fl_Check_Button *SourceIsolines;

/// Widget for target file name
Fl_Output *filename_target;

/// Widget for source file name
Fl_Output *filename_source;

/// Widget for coordinate controls (flip X)
Fl_Check_Button *FlipX;

/// Widget for coordinate controls (flip Y)
Fl_Check_Button *FlipY;

/// Widget for coordinate controls (flip Z)
Fl_Check_Button *FlipZ;

/// Widget for coordinate controls
Fl_Check_Button *DisplayNeurological;

/// Widget for coordinate controls
Fl_Check_Button *DisplayRadiological;

/// Widget for coordinate controls
Fl_Check_Button *DisplayNative;

/// Widget for target viewing mode
Fl_Button *viewTarget;

/// Widget for source viewing mode
Fl_Button *viewSource;

/// Widget for horizontal shutter viewing mode
Fl_Button *viewHShutter;

/// Widget for vertical shutter viewing mode
Fl_Button *viewVShutter;

/// Widget for subtraction viewing mode
Fl_Button *viewSubtraction;

/// Widget for checkerboard viewing mode
Fl_Button *viewCheckerboard;

/// Widget for A over B viewing mode
Fl_Button *viewAoverB;

/// Widget for B over A viewing mode
Fl_Button *viewBoverA;

/// Widget for display viewing mix in shutter viewing mode
Fl_Value_Slider *viewMix;

// Callbacks for image display
static void cb_targetMinMax(Fl_Value_Slider*, void*);
static void cb_sourceMinMax(Fl_Value_Slider*, void*);
static void cb_targetFrame(Fl_Value_Slider*, void*);
static void cb_sourceFrame(Fl_Value_Slider*, void*);
static void cb_lineThickness(Fl_Value_Slider*, void*);
static void cb_speed(Fl_Value_Slider*, void*);
static void cb_TargetIsolines(Fl_Check_Button*, void*);
static void cb_SourceIsolines(Fl_Check_Button*, void*);
static void cb_flipCoordinates(Fl_Check_Button*, void*);
static void cb_displayCoordinates(Fl_Check_Button*, void*);
static void cb_TargetInterpolation(Fl_Menu_*, void*);
static void cb_SourceInterpolation(Fl_Menu_*, void*);
static void cb_TargetColor(Fl_Menu_*, void*);
static void cb_SourceColor(Fl_Menu_*, void*);
static void cb_viewMix(Fl_Value_Slider*, void*);
static void cb_viewTarget(Fl_Button*, void*);
static void cb_viewSource(Fl_Button*, void*);
static void cb_viewHShutter(Fl_Button*, void*);
static void cb_viewVShutter(Fl_Button*, void*);
static void cb_viewSubtraction(Fl_Button*, void*);
static void cb_viewCheckerboard(Fl_Button*, void*);
static void cb_viewAoverB(Fl_Button*, void*);
static void cb_viewBoverA(Fl_Button*, void*);

// Callbacks for image I/O
static void cb_loadTarget(Fl_Button*, void*);
static void cb_loadSource(Fl_Button*, void*);
static void cb_saveSource(Fl_Button*, void*);
static void cb_saveTarget(Fl_Button*, void*);
static void cb_saveScreen(Fl_Button*, void*);

// Callbacks for animation
static void cb_movieStart(Fl_Button*, void*);
static void cb_movieStop(Fl_Button*, void*);
static void cb_playback(void*);
static void cb_savePlayback(Fl_Button*, void *);

public:

/// Add target image
void AddTarget(char *);

/// Add source image
void AddSource(char *);

/// Show the image control window
void ShowImageControlWindow();

/// Update the image control window
void UpdateImageControlWindow();

/// Initialize the image control window
void InitializeImageControlWindow();
