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

#include <sys/types.h>
#include <sys/stat.h>

// Local includes
#include <Fl_RViewUI.h>

#include <bitmaps/target.xbm>
#include <bitmaps/source.xbm>
#include <bitmaps/AoverB.xbm>
#include <bitmaps/BoverA.xbm>
#include <bitmaps/hshutter.xbm>
#include <bitmaps/vshutter.xbm>
#include <bitmaps/subtraction.xbm>
#include <bitmaps/checkerboard.xbm>

#ifdef HAS_VISUALISATION_PANEL

extern Fl_RViewUI  *rviewUI;
extern Fl_RView    *viewer;
extern irtkRView    *rview;

char interpStrings[5][255] = {"nn", "linear", "c1spline", "bspline", "sinc"};

Fl_Menu_Item Fl_RViewUI::menu_targetInterpolationMode[] = {
  {"NN Interpolation",        0, (Fl_Callback*)cb_TargetInterpolation, interpStrings[0], 0, 0, 0, 0, 0},
  {"Linear Interpolation",    0, (Fl_Callback*)cb_TargetInterpolation, interpStrings[1], 0, 0, 0, 0, 0},
  {"C1-spline Interpolation", 0, (Fl_Callback*)cb_TargetInterpolation, interpStrings[2], 0, 0, 0, 0, 0},
  {"B-spline Interpolation",  0, (Fl_Callback*)cb_TargetInterpolation, interpStrings[3], 0, 0, 0, 0, 0},
  {"Sinc Interpolation",      0, (Fl_Callback*)cb_TargetInterpolation, interpStrings[4], 0, 0, 0, 0, 0},
  { 0, 0, 0, 0, 0, 0, 0, 0, 0 }
};

Fl_Menu_Item Fl_RViewUI::menu_sourceInterpolationMode[] = {
  {"NN Interpolation",        0, (Fl_Callback*)cb_SourceInterpolation, interpStrings[0], 0, 0, 0, 0, 0},
  {"Linear Interpolation",    0, (Fl_Callback*)cb_SourceInterpolation, interpStrings[1], 0, 0, 0, 0, 0},
  {"C1-spline Interpolation", 0, (Fl_Callback*)cb_SourceInterpolation, interpStrings[2], 0, 0, 0, 0, 0},
  {"B-spline Interpolation",  0, (Fl_Callback*)cb_SourceInterpolation, interpStrings[3], 0, 0, 0, 0, 0},
  {"Sinc Interpolation",      0, (Fl_Callback*)cb_SourceInterpolation, interpStrings[4], 0, 0, 0, 0, 0},
  { 0, 0, 0, 0, 0, 0, 0, 0, 0 }
};

char colorStrings[11][255] = {"red", "green", "blue", "grey", "inverse", "jacobian", "jacobian expansion", "jacobian contraction", "hotmetal", "rainbow", "custom"};

Fl_Menu_Item Fl_RViewUI::menu_targetColorMode[] = {
  { "Red",      0, (Fl_Callback *)cb_TargetColor, colorStrings[0], 0, 0, 0, 0, 0},
  { "Green",    0, (Fl_Callback *)cb_TargetColor, colorStrings[1], 0, 0, 0, 0, 0},
  { "Blue",     0, (Fl_Callback *)cb_TargetColor, colorStrings[2], 0, 0, 0, 0, 0},
  { "Grey",     0, (Fl_Callback *)cb_TargetColor, colorStrings[3], 0, 0, 0, 0, 0},
  { "Inverse",  0, (Fl_Callback *)cb_TargetColor, colorStrings[4], 0, 0, 0, 0, 0},
  { "Jacobian", 0, (Fl_Callback *)cb_TargetColor, colorStrings[5], 0, 0, 0, 0, 0},
  { "Jacobian / Expansion", 0, (Fl_Callback *)cb_TargetColor, colorStrings[6], 0, 0, 0, 0, 0},
  { "Jacobian / Contracion", 0, (Fl_Callback *)cb_TargetColor, colorStrings[7], 0, 0, 0, 0, 0},
  { "Hotmetal", 0, (Fl_Callback *)cb_TargetColor, colorStrings[8], 0, 0, 0, 0, 0},
  { "Rainbow",  0, (Fl_Callback *)cb_TargetColor, colorStrings[9], 0, 0, 0, 0, 0},
  { "Custom",   0, (Fl_Callback *)cb_TargetColor, colorStrings[10], 0, 0, 0, 0, 0},
  { 0, 0, 0, 0, 0, 0, 0, 0, 0 }
};

Fl_Menu_Item Fl_RViewUI::menu_sourceColorMode[] = {
  { "Red",      0, (Fl_Callback *)cb_SourceColor, colorStrings[0], 0, 0, 0, 0, 0},
  { "Green",    0, (Fl_Callback *)cb_SourceColor, colorStrings[1], 0, 0, 0, 0, 0},
  { "Blue",     0, (Fl_Callback *)cb_SourceColor, colorStrings[2], 0, 0, 0, 0, 0},
  { "Grey",     0, (Fl_Callback *)cb_SourceColor, colorStrings[3], 0, 0, 0, 0, 0},
  { "Inverse",  0, (Fl_Callback *)cb_SourceColor, colorStrings[4], 0, 0, 0, 0, 0},
  { "Jacobian", 0, (Fl_Callback *)cb_SourceColor, colorStrings[5], 0, 0, 0, 0, 0},
  { "Jacobian / Expansion", 0, (Fl_Callback *)cb_SourceColor, colorStrings[6], 0, 0, 0, 0, 0},
  { "Jacobian / Contraction", 0, (Fl_Callback *)cb_SourceColor, colorStrings[7], 0, 0, 0, 0, 0},
  { "Hotmetal", 0, (Fl_Callback *)cb_SourceColor, colorStrings[8], 0, 0, 0, 0, 0},
  { "Rainbow",  0, (Fl_Callback *)cb_SourceColor, colorStrings[9], 0, 0, 0, 0, 0},
  { "Custom",   0, (Fl_Callback *)cb_SourceColor, colorStrings[10], 0, 0, 0, 0, 0},
  { 0, 0, 0, 0, 0, 0, 0, 0, 0 }
};

void Fl_RViewUI::AddTarget(char *filename)
{
  rviewUI->filename_target->value(filename);
}

void Fl_RViewUI::AddSource(char *filename)
{
  rviewUI->filename_source->value(filename);
}

void Fl_RViewUI::cb_loadTarget(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Load target image", "*.{gipl,gipl.Z,hdr,hdr.gz,nii,nii.gz}", "");
  if (filename != NULL) {
    rview->ReadTarget(filename);
    rview->Update();
    viewer->redraw();
    rviewUI->update();
    rviewUI->filename_target->value(filename);
  }
}

void Fl_RViewUI::cb_loadSource(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Load source image", "*.{gipl,gipl.Z,hdr,hdr.gz,nii,nii.gz}", "");
  if (filename != NULL) {
    rview->ReadSource(filename);
    rview->Update();
    viewer->redraw();
    rviewUI->update();
    rviewUI->filename_source->value(filename);
  }
}

void Fl_RViewUI::cb_saveTarget(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Save target image", "*.{gipl,hdr,vtk,nii}", "");
  if (filename != NULL) {
    struct stat buffer;
    if (stat(filename, &buffer) == 0) {
      if (fl_choice("File already exists. Do you want to overwrite the file?", NULL, "No", "Yes") == 2) {
        rview->WriteTarget(filename);
      }
    } else {
      rview->WriteTarget(filename);
    }
  }
}

void Fl_RViewUI::cb_saveSource(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Save source image", "*.{gipl,hdr,vtk,nii}", "");
  if (filename != NULL) {
    struct stat buffer;
    if (stat(filename, &buffer) == 0) {
      if (fl_choice("File already exists. Do you want to overwrite the file?", NULL, "No", "Yes") == 2) {
        rview->WriteSource(filename);
      }
    } else {
      rview->WriteSource(filename);
    }
  }
}

void Fl_RViewUI::cb_movieStart(Fl_Button *, void *)
{
  if (!Fl::has_idle(cb_playback)) Fl::add_idle(cb_playback);
}

void Fl_RViewUI::cb_movieStop(Fl_Button *, void *)
{
  if (Fl::has_idle(cb_playback)) Fl::remove_idle(cb_playback);
}

void Fl_RViewUI::cb_playback(void *)
{
  int t;

  t = rview->GetTargetFrame()+1;
  if (t < 0) t = rview->GetTarget()->GetT()-1;
  if (t >= rview->GetTarget()->GetT()) t = 0;
  rview->SetTargetFrame(round(t));
  rview->SetSourceFrame(round(t));
  rview->Update();
  viewer->redraw();
  rviewUI->update();
}

void Fl_RViewUI::cb_savePlayback(Fl_Button *, void *)
{
  int i, t, s;
  char buffer2[256];

  if (fl_choice("Do you want to save the movie to disk?", NULL, "No", "Yes") == 2) {

    t = rview->GetTargetFrame();
    s = rview->GetSourceFrame();
    for (i = 0; i < rview->GetTarget()->GetT(); i++) {
      rview->SetTargetFrame(i);
      rview->SetSourceFrame(i);

      // Update
      rview->Update();
      sprintf(buffer2, "movie_%.5d.png", i);
      rview->DrawOffscreen(buffer2);

      // Draw on screen to
      viewer->redraw();

      // Force drawing
      Fl::wait(0);

    }
    rview->SetTargetFrame(t);
    rview->SetSourceFrame(s);

  } else {

    t = rview->GetTargetFrame();
    s = rview->GetSourceFrame();
    for (i = 0; i < rview->GetTarget()->GetT(); i++) {
      rview->SetTargetFrame(i);
      rview->SetSourceFrame(i);

      // Update
      rview->Update();
      viewer->redraw();

      // Force drawing
      Fl::wait(0);
    }
    rview->SetTargetFrame(t);
    rview->SetSourceFrame(s);

  }

  // Update
  rview->Update();
  viewer->redraw();

}

void Fl_RViewUI::cb_TargetIsolines(Fl_Check_Button* o, void*)
{
  if (o->value() == 0) rview->DisplayTargetContoursOff();;
  if (o->value() == 1) rview->DisplayTargetContoursOn();;
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_lineThickness(Fl_Value_Slider* o, void*)
{
  rview->SetLineThickness(o->value());
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_SourceIsolines(Fl_Check_Button* o, void*)
{
  if (o->value() == 0) rview->DisplaySourceContoursOff();
  if (o->value() == 1) rview->DisplaySourceContoursOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_TargetInterpolation(Fl_Menu_*, void* v)
{
  if (strcmp((char *)v, "nn") == 0) {
    rview->SetTargetInterpolationMode(Interpolation_NN);
  }
  if (strcmp((char *)v, "linear") == 0) {
    rview->SetTargetInterpolationMode(Interpolation_Linear);
  }
  if (strcmp((char *)v, "c1spline") == 0) {
    rview->SetTargetInterpolationMode(Interpolation_CSpline);
  }
  if (strcmp((char *)v, "bspline") == 0) {
    rview->SetTargetInterpolationMode(Interpolation_BSpline);
  }
  if (strcmp((char *)v, "sinc") == 0) {
    rview->SetTargetInterpolationMode(Interpolation_Sinc);
  }
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_SourceInterpolation(Fl_Menu_*, void* v)
{
  if (strcmp((char *)v, "nn") == 0) {
    rview->SetSourceInterpolationMode(Interpolation_NN);
  }
  if (strcmp((char *)v, "linear") == 0) {
    rview->SetSourceInterpolationMode(Interpolation_Linear);
  }
  if (strcmp((char *)v, "c1spline") == 0) {
    rview->SetSourceInterpolationMode(Interpolation_CSpline);
  }
  if (strcmp((char *)v, "bspline") == 0) {
    rview->SetSourceInterpolationMode(Interpolation_BSpline);
  }
  if (strcmp((char *)v, "sinc") == 0) {
    rview->SetSourceInterpolationMode(Interpolation_Sinc);
  }
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_TargetColor(Fl_Menu_ *, void *v)
{
  if (strcmp((char *)v, "red") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToRed();
  }
  if (strcmp((char *)v, "green") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToGreen();
  }
  if (strcmp((char *)v, "blue") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToBlue();
  }
  if (strcmp((char *)v, "grey") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToLuminance();
  }
  if (strcmp((char *)v, "inverse") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToInverseLuminance();
  }
  if (strcmp((char *)v, "jacobian") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToJacobian();
  }
  if (strcmp((char *)v, "jacobian expansion") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToJacobianExpansion();
  }
  if (strcmp((char *)v, "jacobian contraction") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToJacobianContraction();
  }
  if (strcmp((char *)v, "hotmetal") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToHotMetal();
  }
  if (strcmp((char *)v, "rainbow") == 0) {
    rview->GetTargetLookupTable()->SetColorModeToRainbow();
  }
  if (strcmp((char *)v, "custom") == 0) {
    char *filename = fl_file_chooser("Load lookup table", "*.lut", "");
    if (filename != NULL) rview->GetTargetLookupTable()->Read(filename);
  }
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_SourceColor(Fl_Menu_ *, void *v)
{
  if (strcmp((char *)v, "red") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToRed();
  }
  if (strcmp((char *)v, "green") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToGreen();
  }
  if (strcmp((char *)v, "blue") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToBlue();
  }
  if (strcmp((char *)v, "grey") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToLuminance();
  }
  if (strcmp((char *)v, "inverse") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToInverseLuminance();
  }
  if (strcmp((char *)v, "jacobian") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToJacobian();
  }
  if (strcmp((char *)v, "jacobian expansion") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToJacobianExpansion();
  }
  if (strcmp((char *)v, "jacobian contraction") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToJacobianContraction();
  }
  if (strcmp((char *)v, "hotmetal") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToHotMetal();
  }
  if (strcmp((char *)v, "rainbow") == 0) {
    rview->GetSourceLookupTable()->SetColorModeToRainbow();
  }
  if (strcmp((char *)v, "custom") == 0) {
    char *filename = fl_file_chooser("Load lookup table", "*.lut", "");
    if (filename != NULL) rview->GetSourceLookupTable()->Read(filename);
  }
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_targetMinMax(Fl_Value_Slider*, void*)
{
  if (rview->GetViewMode() != View_Subtraction) {
    rview->SetDisplayMinTarget(rviewUI->targetMin->value());
    rview->SetDisplayMaxTarget(rviewUI->targetMax->value());
  } else {
    rview->SetDisplayMinSubtraction(rviewUI->targetMin->value());
    rview->SetDisplayMaxSubtraction(rviewUI->targetMax->value());
  }
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_sourceMinMax(Fl_Value_Slider*, void*)
{
  rview->SetDisplayMinSource(rviewUI->sourceMin->value());
  rview->SetDisplayMaxSource(rviewUI->sourceMax->value());
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_targetFrame(Fl_Value_Slider*, void*)
{
  rview->SetTargetFrame(round(rviewUI->targetFrame->value()));
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_sourceFrame(Fl_Value_Slider*, void*)
{
  rview->SetSourceFrame(round(rviewUI->sourceFrame->value()));
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_flipCoordinates(Fl_Check_Button*, void*)
{
  if (rviewUI->FlipX->value() == 0) rview->FlipXOff();
  if (rviewUI->FlipX->value() == 1) rview->FlipXOn();
  if (rviewUI->FlipY->value() == 0) rview->FlipYOff();
  if (rviewUI->FlipY->value() == 1) rview->FlipYOn();
  if (rviewUI->FlipZ->value() == 0) rview->FlipZOff();
  if (rviewUI->FlipZ->value() == 1) rview->FlipZOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_displayCoordinates(Fl_Check_Button*, void*)
{
  if (rviewUI->DisplayNeurological->value() == 1) rview->SetDisplayMode(Neurological);
  if (rviewUI->DisplayRadiological->value() == 1) rview->SetDisplayMode(Radiological);
  if (rviewUI->DisplayNative->value() == 1) rview->SetDisplayMode(Native);
  rview->Reset();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewMix(Fl_Value_Slider *o, void *)
{
  rview->SetViewMix(o->value());
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewTarget(Fl_Button*, void*)
{
  if (rview->GetViewMode() == View_Subtraction) {
    rviewUI->targetMin->minimum(rview->GetTargetMin());
    rviewUI->targetMin->maximum(rview->GetTargetMax());
    rviewUI->targetMax->minimum(rview->GetTargetMin());
    rviewUI->targetMax->maximum(rview->GetTargetMax());
    rviewUI->targetMin->value(rview->GetDisplayMinTarget());
    rviewUI->targetMax->value(rview->GetDisplayMaxTarget());
  }
  rview->SetViewMode(View_A);
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewSource(Fl_Button*, void*)
{
  if (rview->GetViewMode() == View_Subtraction) {
    rviewUI->targetMin->minimum(rview->GetTargetMin());
    rviewUI->targetMin->maximum(rview->GetTargetMax());
    rviewUI->targetMax->minimum(rview->GetTargetMin());
    rviewUI->targetMax->maximum(rview->GetTargetMax());
    rviewUI->targetMin->value(rview->GetDisplayMinTarget());
    rviewUI->targetMax->value(rview->GetDisplayMaxTarget());
  }
  rview->SetViewMode(View_B);
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewHShutter(Fl_Button*, void*)
{
  if (rview->GetViewMode() == View_Subtraction) {
    rviewUI->targetMin->minimum(rview->GetTargetMin());
    rviewUI->targetMin->maximum(rview->GetTargetMax());
    rviewUI->targetMax->minimum(rview->GetTargetMin());
    rviewUI->targetMax->maximum(rview->GetTargetMax());
    rviewUI->targetMin->value(rview->GetDisplayMinTarget());
    rviewUI->targetMax->value(rview->GetDisplayMaxTarget());
  }
  rview->SetViewMode(View_HShutter);
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewVShutter(Fl_Button*, void*)
{
  if (rview->GetViewMode() == View_Subtraction) {
    rviewUI->targetMin->minimum(rview->GetTargetMin());
    rviewUI->targetMin->maximum(rview->GetTargetMax());
    rviewUI->targetMax->minimum(rview->GetTargetMin());
    rviewUI->targetMax->maximum(rview->GetTargetMax());
    rviewUI->targetMin->value(rview->GetDisplayMinTarget());
    rviewUI->targetMax->value(rview->GetDisplayMaxTarget());
  }
  rview->SetViewMode(View_VShutter);
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewSubtraction(Fl_Button*, void*)
{
  if (rview->GetViewMode() != View_Subtraction) {
    rviewUI->targetMin->minimum(rview->GetSubtractionMin());
    rviewUI->targetMin->maximum(rview->GetSubtractionMax());
    rviewUI->targetMax->minimum(rview->GetSubtractionMin());
    rviewUI->targetMax->maximum(rview->GetSubtractionMax());
    rviewUI->targetMin->value(rview->GetDisplayMinSubtraction());
    rviewUI->targetMax->value(rview->GetDisplayMaxSubtraction());
  }
  rview->SetViewMode(View_Subtraction);
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewCheckerboard(Fl_Button *, void*)
{
  if (rview->GetViewMode() == View_Subtraction) {
    rviewUI->targetMin->minimum(rview->GetTargetMin());
    rviewUI->targetMin->maximum(rview->GetTargetMax());
    rviewUI->targetMax->minimum(rview->GetTargetMin());
    rviewUI->targetMax->maximum(rview->GetTargetMax());
    rviewUI->targetMin->value(rview->GetDisplayMinTarget());
    rviewUI->targetMax->value(rview->GetDisplayMaxTarget());
  }
  rview->SetViewMode(View_Checkerboard);
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewAoverB(Fl_Button *, void*)
{
  if (rview->GetViewMode() == View_Subtraction) {
    rviewUI->targetMin->minimum(rview->GetTargetMin());
    rviewUI->targetMin->maximum(rview->GetTargetMax());
    rviewUI->targetMax->minimum(rview->GetTargetMin());
    rviewUI->targetMax->maximum(rview->GetTargetMax());
    rviewUI->targetMin->value(rview->GetDisplayMinTarget());
    rviewUI->targetMax->value(rview->GetDisplayMaxTarget());
  }
  rview->SetViewMode(View_AoverB);
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewBoverA(Fl_Button *, void*)
{
  if (rview->GetViewMode() == View_Subtraction) {
    rviewUI->targetMin->minimum(rview->GetTargetMin());
    rviewUI->targetMin->maximum(rview->GetTargetMax());
    rviewUI->targetMax->minimum(rview->GetTargetMin());
    rviewUI->targetMax->maximum(rview->GetTargetMax());
    rviewUI->targetMin->value(rview->GetDisplayMinTarget());
    rviewUI->targetMax->value(rview->GetDisplayMaxTarget());
  }
  rview->SetViewMode(View_BoverA);
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::ShowImageControlWindow()
{
  targetMin->minimum(rview->GetTargetMin());
  targetMin->maximum(rview->GetTargetMax());
  targetMax->minimum(rview->GetTargetMin());
  targetMax->maximum(rview->GetTargetMax());
  targetMin->value(rview->GetDisplayMinSource());
  targetMax->value(rview->GetDisplayMaxSource());

  // Set up correct min and max values for all source widgets
  sourceMin->minimum(rview->GetSourceMin());
  sourceMin->maximum(rview->GetSourceMax());
  sourceMax->minimum(rview->GetSourceMin());
  sourceMax->maximum(rview->GetSourceMax());
  sourceMin->value(rview->GetDisplayMinSource());
  sourceMax->value(rview->GetDisplayMaxSource());

}

void Fl_RViewUI::UpdateImageControlWindow()
{
  // Update interpolation menus
  switch (rview->GetTargetInterpolationMode()) {
  case Interpolation_NN:
    targetInterpolationMode->value(0);
    break;
  case Interpolation_Linear:
    targetInterpolationMode->value(1);
    break;
  case Interpolation_CSpline:
    targetInterpolationMode->value(2);
    break;
  case Interpolation_BSpline:
    targetInterpolationMode->value(3);
    break;
  case Interpolation_Sinc:
    targetInterpolationMode->value(4);
    break;
  default:
    break;
  }
  switch (rview->GetSourceInterpolationMode()) {
  case Interpolation_NN:
    sourceInterpolationMode->value(0);
    break;
  case Interpolation_Linear:
    sourceInterpolationMode->value(1);
    break;
  case Interpolation_CSpline:
    sourceInterpolationMode->value(2);
    break;
  case Interpolation_BSpline:
    sourceInterpolationMode->value(3);
    break;
  case Interpolation_Sinc:
    sourceInterpolationMode->value(4);
    break;
  default:
    break;
  }
  switch (rview->GetTargetLookupTable()->GetColorMode()) {
  case ColorMode_Red:
    targetColorMode->value(0);
    break;
  case ColorMode_Green:
    targetColorMode->value(1);
    break;
  case ColorMode_Blue:
    targetColorMode->value(2);
    break;
  case ColorMode_Luminance:
    targetColorMode->value(3);
    break;
  case ColorMode_InverseLuminance:
    targetColorMode->value(4);
    break;
  case ColorMode_Jacobian:
    targetColorMode->value(5);
    break;
  case ColorMode_JacobianExpansion:
    targetColorMode->value(6);
    break;
  case ColorMode_JacobianContraction:
    targetColorMode->value(7);
    break;
  case ColorMode_HotMetal:
    targetColorMode->value(8);
    break;
  case ColorMode_Rainbow:
    targetColorMode->value(9);
    break;
  case ColorMode_Custom:
    targetColorMode->value(10);
    break;
  default:
    break;
  }
  switch (rview->GetSourceLookupTable()->GetColorMode()) {
  case ColorMode_Red:
    sourceColorMode->value(0);
    break;
  case ColorMode_Green:
    sourceColorMode->value(1);
    break;
  case ColorMode_Blue:
    sourceColorMode->value(2);
    break;
  case ColorMode_Luminance:
    sourceColorMode->value(3);
    break;
  case ColorMode_InverseLuminance:
    sourceColorMode->value(4);
    break;
  case ColorMode_Jacobian:
    sourceColorMode->value(5);
    break;
  case ColorMode_JacobianExpansion:
    sourceColorMode->value(6);
    break;
  case ColorMode_JacobianContraction:
    sourceColorMode->value(7);
    break;
  case ColorMode_HotMetal:
    sourceColorMode->value(8);
    break;
  case ColorMode_Rainbow:
    sourceColorMode->value(9);
    break;
  case ColorMode_Custom:
    sourceColorMode->value(10);
    break;
  default:
    break;
  }

  // Update isoline menus
  rviewUI->TargetIsolines->value(rview->GetDisplayTargetContours());
  rviewUI->SourceIsolines->value(rview->GetDisplaySourceContours());

  // Update Line Thickness
  rviewUI->lineThickness->value(rview->GetLineThickness());

  // Update viewing mode
  rviewUI->DisplayNeurological->value(0);
  rviewUI->DisplayRadiological->value(0);
  rviewUI->DisplayNative->value(0);

  switch (rview->GetDisplayMode()) {
  case Neurological:
    rviewUI->DisplayNeurological->value(1);
    break;
  case Radiological:
    rviewUI->DisplayRadiological->value(1);
    break;
  case Native:
    rviewUI->DisplayNative->value(1);
    break;
  }

  // Update viewing mode
  rviewUI->viewTarget->value(0);
  rviewUI->viewSource->value(0);
  rviewUI->viewHShutter->value(0);
  rviewUI->viewVShutter->value(0);
  rviewUI->viewSubtraction->value(0);
  rviewUI->viewCheckerboard->value(0);

  switch (rview->GetViewMode()) {
  case View_A:
    rviewUI->viewTarget->value(1);
    break;
  case View_B:
    rviewUI->viewSource->value(1);
    break;
  case View_Checkerboard:
    rviewUI->viewCheckerboard->value(1);
    break;
  case View_Subtraction:
    rviewUI->viewSubtraction->value(1);
    break;
  case View_HShutter:
    rviewUI->viewHShutter->value(1);
    break;
  case View_VShutter:
    rviewUI->viewVShutter->value(1);
    break;
  case View_AoverB:
    rviewUI->viewAoverB->value(1);
    break;
  case View_BoverA:
    rviewUI->viewBoverA->value(1);
    break;
  }

  // Update viewing mix
  rviewUI->viewMix->value(rview->GetViewMix());

  // Set up correct min and max values for all target widgets
  if (rview->GetViewMode() != View_Subtraction) {
    targetMin->minimum(rview->GetTargetMin());
    targetMin->maximum(rview->GetTargetMax());
    targetMax->minimum(rview->GetTargetMin());
    targetMax->maximum(rview->GetTargetMax());
    targetMin->value(rview->GetDisplayMinTarget());
    targetMax->value(rview->GetDisplayMaxTarget());
  } else {
    targetMin->minimum(rview->GetSubtractionMin());
    targetMin->maximum(rview->GetSubtractionMax());
    targetMax->minimum(rview->GetSubtractionMin());
    targetMax->maximum(rview->GetSubtractionMax());
    targetMin->value(rview->GetDisplayMinSubtraction());
    targetMax->value(rview->GetDisplayMaxSubtraction());
  }

  // Set up correct min and max values for all source widgets
  sourceMin->minimum(rview->GetSourceMin());
  sourceMin->maximum(rview->GetSourceMax());
  sourceMax->minimum(rview->GetSourceMin());
  sourceMax->maximum(rview->GetSourceMax());
  sourceMin->value(rview->GetDisplayMinSource());
  sourceMax->value(rview->GetDisplayMaxSource());

  // Set up correct min and max values for time frames
  targetFrame->minimum(0);
  targetFrame->maximum(rview->GetTarget()->GetT()-1);
  targetFrame->value(rview->GetTargetFrame());
  sourceFrame->minimum(0);
  sourceFrame->maximum(rview->GetSource()->GetT()-1);
  sourceFrame->value(rview->GetSourceFrame());
}

void Fl_RViewUI::InitializeImageControlWindow()
{
  {
    // Create target controls
    Fl_Group* o = new Fl_Group(0, 30, 400, 200, "Target");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Value_Slider2* o = targetMin =
                             new Fl_Value_Slider2(10, 50, 380, 20, "Min. greyvalue:");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->callback((Fl_Callback*)cb_targetMinMax);
      o->align(133);
    }
    {
      Fl_Value_Slider2* o = targetMax =
                             new Fl_Value_Slider2(10, 90, 380, 20, "Max. greyvalue:");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->callback((Fl_Callback*)cb_targetMinMax);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    {
      Fl_Value_Slider* o = targetFrame = new Fl_Value_Slider(10, 130, 380, 20, "Frame no:");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->step(1);
      o->callback((Fl_Callback*)cb_targetFrame);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    {
      Fl_Choice* o = targetInterpolationMode = new Fl_Choice(10, 160, 140, 30);
      o->align(FL_ALIGN_BOTTOM);
      o->menu(menu_targetInterpolationMode);
    }
    {
      Fl_Choice* o = targetColorMode = new Fl_Choice(170, 160, 120, 30);
      o->align(FL_ALIGN_BOTTOM);
      o->menu(menu_targetColorMode);
    }
    {
      Fl_Check_Button* o = TargetIsolines = new Fl_Check_Button(310, 160, 80, 30, "Isolines");
      o->callback((Fl_Callback*)cb_TargetIsolines);
    }
    {
      Fl_Output* o = filename_target = new Fl_Output(70, 200, 320, 20, "Image = ");
      o->box(FL_FLAT_BOX);
    }

    o->end(); // End of target controls
  }
  {
    // Create source controls
    Fl_Group* o = new Fl_Group(0, 230, 400, 200, "Source");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Value_Slider2* o = sourceMin = new Fl_Value_Slider2(10, 250, 380, 20, "Min. greyvalue:");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->callback((Fl_Callback*)cb_sourceMinMax);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    {
      Fl_Value_Slider2* o = sourceMax = new Fl_Value_Slider2(10, 290, 380, 20, "Max. greyvalue:");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->callback((Fl_Callback*)cb_sourceMinMax);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    {
      Fl_Value_Slider* o = sourceFrame = new Fl_Value_Slider(10, 330, 380, 20, "Frame no:");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->step(1);
      o->callback((Fl_Callback*)cb_sourceFrame);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    {
      Fl_Choice* o = sourceInterpolationMode = new Fl_Choice(10, 360, 140, 30);
      o->align(FL_ALIGN_BOTTOM);
      o->menu(menu_sourceInterpolationMode);
    }
    {
      Fl_Choice* o = sourceColorMode = new Fl_Choice(170, 360, 120, 30);
      o->align(FL_ALIGN_BOTTOM);
      o->menu(menu_sourceColorMode);
    }
    {
      Fl_Check_Button* o = SourceIsolines = new Fl_Check_Button(310, 360, 80, 30, "Isolines");
      o->callback((Fl_Callback*)cb_SourceIsolines);
    }
    {
      Fl_Output* o = filename_source = new Fl_Output(70, 400, 320, 20, "Image = ");
      o->box(FL_FLAT_BOX);
    }
    o->end(); // End of source controls
  }

  {
    // Create coordinate controls
    Fl_Group* o = new Fl_Group(0, 430, 400, 40, "Coordinates");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Check_Button* o = DisplayNeurological =
                             new Fl_Check_Button(25, 435, 100, 30, "Neurological");
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->type(FL_RADIO_BUTTON);
      o->callback((Fl_Callback*)cb_displayCoordinates);
    }
    {
      Fl_Check_Button* o = DisplayRadiological =
                             new Fl_Check_Button(150, 435, 100, 30, "Radiological");
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->type(FL_RADIO_BUTTON);
      o->callback((Fl_Callback*)cb_displayCoordinates);
    }
    {
      Fl_Check_Button* o = DisplayNative =
                             new Fl_Check_Button(275, 435, 100, 30, "Native");
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->type(FL_RADIO_BUTTON);
      o->callback((Fl_Callback*)cb_displayCoordinates);
    }
    o->end(); // End of coordinate controls
  }

  {
    // Create display controls
    Fl_Group* o = new Fl_Group(0, 470, 400, 190, "Display");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Button* o = viewTarget = new Fl_Button(40, 490, 50, 50);
      o->type(102);
      o->value(1);
      o->callback((Fl_Callback*)cb_viewTarget);
      Fl_Bitmap *p = new Fl_Bitmap(target_bits, 50, 50);
      p->label(o);
    }
    {
      Fl_Button* o = viewSource = new Fl_Button(40, 550, 50, 50);
      o->type(102);
      o->callback((Fl_Callback*)cb_viewSource);
      Fl_Bitmap *p = new Fl_Bitmap(source_bits, 50, 50);
      p->label(o);
    }
    {
      Fl_Button* o = viewHShutter = new Fl_Button(130, 490, 50, 50);
      o->type(102);
      o->callback((Fl_Callback*)cb_viewHShutter);
      Fl_Bitmap *p = new Fl_Bitmap(hshutter_bits, 50, 50);
      p->label(o);
    }
    {
      Fl_Button* o = viewVShutter = new Fl_Button(130, 550, 50, 50);
      o->type(102);
      o->callback((Fl_Callback*)cb_viewVShutter);
      Fl_Bitmap *p = new Fl_Bitmap(vshutter_bits, 50, 50);
      p->label(o);
    }
    {
      Fl_Button* o = viewSubtraction = new Fl_Button(220, 490, 50, 50);
      o->type(102);
      o->callback((Fl_Callback*)cb_viewSubtraction);
      Fl_Bitmap *p = new Fl_Bitmap(subtraction_bits, 50, 50);
      p->label(o);
    }
    {
      Fl_Button* o = viewCheckerboard = new Fl_Button(220, 550, 50, 50);
      o->type(102);
      o->callback((Fl_Callback*)cb_viewCheckerboard);
      Fl_Bitmap *p = new Fl_Bitmap(checkerboard_bits, 50, 50);
      p->label(o);
    }
    {
      Fl_Button* o = viewAoverB = new Fl_Button(310, 490, 50, 50);
      o->type(102);
      o->callback((Fl_Callback*)cb_viewAoverB);
      Fl_Bitmap *p = new Fl_Bitmap(AoverB_bits, 50, 50);
      p->label(o);
    }
    {
      Fl_Button* o = viewBoverA = new Fl_Button(310, 550, 50, 50);
      o->type(102);
      o->callback((Fl_Callback*)cb_viewBoverA);
      Fl_Bitmap *p = new Fl_Bitmap(BoverA_bits, 50, 50);
      p->label(o);
    }
    {
      Fl_Value_Slider* o = viewMix = new Fl_Value_Slider(10, 630, 380, 20, "Display mix:");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->step(0.01);
      o->value(0.5);
      o->callback((Fl_Callback*)cb_viewMix);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    o->end(); // End of display controls
  }
}

#endif
