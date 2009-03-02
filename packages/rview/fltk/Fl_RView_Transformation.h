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
// Members for transformation display and manipulation
//

// Widget for transformation browsing
Fl_Simple_Browser *transformationBrowser;

// Widget for transformation file name
Fl_Output *info_trans_filename;

// Widget for transformation details
Fl_Hold_Browser *info_trans_details;

/// Widget for applying transformation
Fl_Button *transformApply;

/// Widget for inverting transformation
Fl_Button *transformInvert;

/// Widget for display control of deformation points
Fl_Check_Button *viewDeformationPoints;

/// Widget for display control of deformation grid
Fl_Check_Button *viewDeformationGrid;

/// Widget for display control of deformation arrows
Fl_Check_Button *viewDeformationArrows;

/// Widget for display control of maximum deformation value
Fl_Value_Slider *deformationMax;

/// Widget for display control of minimum deformation value
Fl_Value_Slider *deformationMin;

/// Widget for display control of blending deformation value
Fl_Value_Slider *deformationBlending;

/// Widget for display control of deformation grid resolution
Fl_Value_Slider *viewDeformationGridResolution;

/// Widget to edit transformation
Fl_Window *editTransformation;

/// Widget for transformation parameters
Fl_Valuator **transformationValuator;

// Callbacks for transformations
static void cb_browseTransformation(Fl_Browser*, void*);
static void cb_deleteTransformation(Fl_Button*, void*);
static void cb_applyTransformation(Fl_Button*, void*);
static void cb_invertTransformation(Fl_Button*, void*);
static void cb_moveupTransformation(Fl_Button*, void*);
static void cb_movedownTransformation(Fl_Button*, void*);
static void cb_editTransformation(Fl_Button*, void*);
static void cb_editTransformationApply(Fl_Button*, void*);
static void cb_editTransformationUpdate(Fl_Valuator*, void*);
static void cb_editTransformationReset(Fl_Button*, void*);
static void cb_viewDeformationGridResolution(Fl_Slider*, void*);
static void cb_viewDeformationGrid(Fl_Button*, void*);
static void cb_viewDeformationPoints(Fl_Button*, void*);
static void cb_viewDeformationArrows(Fl_Button*, void*);
static void cb_loadTransformation(Fl_Button*, void*);
static void cb_saveTransformation(Fl_Button*, void*);
static void cb_movieTransformation(Fl_Button*, void*);

public:

/// Add a transformation to the transformation browser
void AddTransformation(char *);

/// Show the transformation control window
void ShowTransformationControlWindow();

/// Update the transformation control window
void UpdateTransformationControlWindow();

/// Initialize the transformation control window
void InitializeTransformationControlWindow();

/// Initialize the transformation editor window
void InitializeTransformationEditor();
