/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _FL_RVIEWUI_H

#define _FL_RVIEWUI_H

#define HAS_VISUALISATION_PANEL
#define HAS_TRANSFORMATION_PANEL
#define HAS_LANDMARK_PANEL
#define HAS_SEGMENTATION_PANEL

// FLTK includes
#include <FL/Fl.H>
#include <FL/Fl_Bitmap.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Browser.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Radio_Button.H>
#include <FL/Fl_Roller.H>
#include <FL/Fl_Choice.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Multi_Browser.H>
#include <FL/Fl_Select_Browser.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/Fl_Hold_Browser.H>
#include <FL/Fl_Check_Browser.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Multiline_Output.H>
#include <FL/Fl_Text_Display.H>
#include <FL/Fl_Return_Button.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Valuator.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Color_Chooser.H>
#include <FL/forms.H>

// FLU
#include <Fl_Simple_Group.h>
#include <Fl_Simple_Browser.h>

// Local includes
#include <Fl_RView.h>
#include <Fl_RView_Histogram.h>

// Slider ranges for affine DOFs
#define _GLOBAL_TRANSLATION_MIN  -80
#define _GLOBAL_TRANSLATION_MAX   80
#define _GLOBAL_SCALING_MIN        1
#define _GLOBAL_SCALING_MAX      200
#define _GLOBAL_ROTATION_MIN    -180
#define _GLOBAL_ROTATION_MAX     180

class Fl_RViewUI
{

  friend class Fl_RView;
  friend class Fl_HistogramWindow;

  //
  // Members for main menu
  //

  /// Widget for main window
  Fl_Window *mainWindow;

  /// Widget for zoom
  Fl_Value_Slider *zoom;

  /// Widget for slice position
  Fl_Roller *sliceX;
  Fl_Roller *sliceY;
  Fl_Roller *sliceZ;

  /// Widget for information
  Fl_Output *info_voxel;
  Fl_Output *info_world;
  Fl_Output *info_target;
  Fl_Output *info_source;
  Fl_Output *info_segmentation;
  Fl_Check_Button *info_snap_to_grid;
  Fl_Check_Button *info_cursor;

  /// Widget for tabs
  Fl_Tabs *tab_menu;

  //
  // Callback methods
  //

  // Callbacks for main menu
  static void cb_SubColor(Fl_Button *, void *);
  static void cb_ViewMode(Fl_Button *, void *);
  static void cb_resliceX(Fl_Roller*, void*);
  static void cb_resliceY(Fl_Roller*, void*);
  static void cb_resliceZ(Fl_Roller*, void*);
  static void cb_zoom(Fl_Value_Slider*, void*);
  static void cb_viewCursor(Fl_Check_Button*, void*);
  static void cb_snapGrid(Fl_Check_Button*, void*);
  static void cb_viewCursorMode(Fl_Button*, void*);

public:

  /// Constructor
  Fl_RViewUI();

  /// Show the windows
  void show();

  /// Update the windows
  void update();

  /// Initialize the main window
  void InitializeMainWindow();

  /// Initialize the cursor control window
  void InitializeCursorControlWindow();

#ifdef HAS_VISUALISATION_PANEL
#include <Fl_RView_Visualisation.h>
#endif

#ifdef HAS_TRANSFORMATION_PANEL
#include <Fl_RView_Transformation.h>
#endif

#ifdef HAS_REGISTRATION_PANEL
#include <Fl_RView_Registration.h>
#endif

#ifdef HAS_LANDMARK_PANEL
#include <Fl_RView_Landmark.h>
#endif

#ifdef HAS_SEGMENTATION_PANEL
#include <Fl_RView_Segmentation.h>
#endif

};

#endif
