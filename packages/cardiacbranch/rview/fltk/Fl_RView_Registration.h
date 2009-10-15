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
// Members for registration display and manipulation
//

Fl_Value_Slider *targetBlurring;
Fl_Value_Slider *sourceBlurring;
Fl_Value_Slider *targetResolution;
Fl_Value_Slider *sourceResolution;
Fl_Value_Slider *targetPadding;
Fl_Counter *numberOfIterations;
Fl_Counter *numberOfBins;
Fl_Counter *numberOfSteps;
Fl_Counter *lengthOfSteps;
Fl_Choice  *similarityMeasure;
Fl_Check_Button *level[5];
Fl_Return_Button *runRegistration;

/// Similarity measure menu
static Fl_Menu_Item menu_similarityMeasure[];

/// Callback for reading paramter
static void cb_loadParameter(Fl_Menu_Bar*, void*);

/// Callback for reading paramter
static void cb_saveParameter(Fl_Menu_Bar*, void*);

/// Callback for parameter updates
static void cb_similarityMeasure(Fl_Widget *, char *) ;

/// Callback for parameter updates
static void cb_numberOfBins(Fl_Counter*, void*);

/// Callback for parameter updates
static void cb_numberOfIterations(Fl_Counter*, void*);

/// Callback for parameter updates
static void cb_numberOfSteps(Fl_Counter*, void*);

/// Callback for parameter updates
static void cb_lengthOfSteps(Fl_Counter*, void*);

// Callbacks for parameter updates
static void cb_numberOfLevels(Fl_Button *, char *);

// Callbacks for parameter updates
static void cb_targetBlurring(Fl_Value_Slider *, void *);

// Callbacks for parameter updates
static void cb_targetResolution(Fl_Value_Slider *, void *);

/// Callback for registration target intensity padding
static void cb_targetPadding(Fl_Value_Slider *, void *);

// Callbacks for parameter updates
static void cb_sourceBlurring(Fl_Value_Slider *, void *);

// Callbacks for parameter updates
static void cb_sourceResolution(Fl_Value_Slider *, void *);

// Callbacks for registrations
static void cb_startRegistration(Fl_Button*, void*);

// Callbacks for registrations
static void cb_stopRegistration(Fl_Button*, void*);

// Callbacks for registrations
static void cb_guessRegistrationParameter(Fl_Button*, void*);

public:

/// Update the registration control window
void UpdateRegistrationControlWindow();

/// Initialize the registration control window
void InitializeRegistrationControlWindow();

