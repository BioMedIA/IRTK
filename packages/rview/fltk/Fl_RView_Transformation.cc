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

#ifdef HAS_TRANSFORMATION_PANEL

// Bitmaps
//#include <bitmaps/grid.xbm>
//#include <bitmaps/arrow.xbm>
//#include <bitmaps/point.xbm>

char deformationPropertyStrings[5][255] = { "none", "displacement", "jacobian", "jacobian_expansion", "jacobian_contraction" };

Fl_Menu_Item Fl_RViewUI::menu_deformationProperty[] = {
  {"None",            				0, (Fl_Callback*)cb_DeformationProperty, deformationPropertyStrings[0]},
  {"Displacement",    				0, (Fl_Callback*)cb_DeformationProperty, deformationPropertyStrings[1]},
  {"Jacobian",        				0, (Fl_Callback*)cb_DeformationProperty, deformationPropertyStrings[2]},
  {"Jacobian (expansion)",    0, (Fl_Callback*)cb_DeformationProperty, deformationPropertyStrings[3]},
  {"Jacobian (contraction)",  0, (Fl_Callback*)cb_DeformationProperty, deformationPropertyStrings[4]},
  {0}
};

extern Fl_RViewUI  *rviewUI;
extern Fl_RView    *viewer;
extern irtkRView    *rview;

void Fl_RViewUI::AddTransformation(char *filename)
{
  int i;

  // Update browser and slider
  rviewUI->transformationBrowser->add(filename);
  rviewUI->transformationBrowser->value(1);

  char *buffer = strdup(rviewUI->transformationBrowser->text(rviewUI->transformationBrowser->value()));
  rview->ReadTransformation(buffer);
  rviewUI->info_trans_filename->value(buffer);
  free(buffer);

  // Update transformation valuator
  irtkTransformation *transform = rview->GetTransformation();
  for (i = 0; i<transform->NumberOfDOFs(); i++) {
    rviewUI->transformationValuator[i]->value(transform->Get(i));
  }

}

void Fl_RViewUI::cb_loadTransformation(Fl_Button *, void *v)
{
  int i;
  char *filename = fl_file_chooser("Load transformation", "*.{dof,dof.gz}", "");

  if (filename != NULL) {
    rview->ReadTransformation(filename);

    // Update browser and info
    rviewUI->transformationBrowser->add(filename);
    rviewUI->transformationBrowser->deselect();
    rviewUI->transformationBrowser->select(rviewUI->transformationBrowser->size());
    rviewUI->info_trans_filename->value(filename);
    // Update
    rview->Update();
    viewer->redraw();
  }

  // Update transformation valuator
  irtkTransformation *transform = rview->GetTransformation();
  for (i = 0; i<transform->NumberOfDOFs(); i++) {
    rviewUI->transformationValuator[i]->value(transform->Get(i));
  }
}

void Fl_RViewUI::cb_saveTransformation(Fl_Button *, void *v)
{
  char *filename = fl_file_chooser("Save transformation", "*.{dof,dof.gz}", "");
  if (filename != NULL) {
    struct stat buffer;
    if (stat(filename, &buffer) == 0) {
      if (fl_choice("File already exists. Do you want to overwrite the file?", NULL, "No", "Yes") == 2) {
        rview->WriteTransformation(filename);
      }
    } else {
      rview->WriteTransformation(filename);
    }
  }
}

void Fl_RViewUI::cb_browseTransformation(Fl_Browser* o, void* v)
{
  if (Fl::event_key() == FL_BackSpace) {
    int n;

    if (rviewUI->transformationBrowser->value() != 0) {
      n = rviewUI->transformationBrowser->value();

      if (rviewUI->transformationBrowser->size() == 1) {
        fl_alert("Can't delete last transformation");
      } else {
        rviewUI->transformationBrowser->remove(n);
        if (n > 1) n--;
        rviewUI->transformationBrowser->select(n);
        char *buffer = strdup(rviewUI->transformationBrowser->text(n));
        rview->ReadTransformation(buffer);
        rviewUI->info_trans_filename->value(buffer);
        free(buffer);
      }
    }
  } else {
    if (o->value() != 0) {
      char *buffer = strdup(o->text(o->value()));
      rview->ReadTransformation(buffer);
      rviewUI->info_trans_filename->value(buffer);
      free(buffer);
    }
  }

  // Update
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_deleteTransformation(Fl_Button* o, void* v)
{
  int n;

  if (rviewUI->transformationBrowser->value() != 0) {
    n = rviewUI->transformationBrowser->value();

    if (rviewUI->transformationBrowser->size() == 1) {
      fl_alert("Can't delete last transformation");
    } else {
      rviewUI->transformationBrowser->remove(n);
      if (n > 1) n--;
      rviewUI->transformationBrowser->select(n);
      char *buffer = strdup(rviewUI->transformationBrowser->text(n));
      rview->ReadTransformation(buffer);
      rviewUI->info_trans_filename->value(buffer);
      free(buffer);
    }
  }
  o->clear();

  // Update
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_editTransformation(Fl_Button* o, void* v)
{
  rviewUI->editTransformation->show();
  o->clear();

  // Update
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}


void Fl_RViewUI::cb_editTransformationApply(Fl_Button* o, void* v)
{
  rviewUI->editTransformation->hide();
}

void Fl_RViewUI::cb_editTransformationUpdate(Fl_Valuator* o, void* v)
{
  // Get transformation
  irtkTransformation *transform = rview->GetTransformation();

  // Set values to selected values
  if ((long)v < transform->NumberOfDOFs()) transform->Put((long)v, o->value());

  // Update
  rview->SourceUpdateOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_editTransformationReset(Fl_Button* o, void* v)
{
  if (fl_choice("Do you really want to reset transformation?", NULL, "No", "Yes") == 2) {
    int i;

    // Get transformation
    irtkTransformation *transform = rview->GetTransformation();

    // Reset values to identity
    for (i=0; i<transform->NumberOfDOFs(); i++) transform->Put(i, 0);
    for (i=6; i<transform->NumberOfDOFs(); i++) transform->Put(i, 100);
    for (i=9; i<transform->NumberOfDOFs(); i++) transform->Put(i, 0);
    for (i=0; i<transform->NumberOfDOFs(); i++) {
      rviewUI->transformationValuator[i]->value(transform->Get(i));
    }
    // Update
    rview->SourceUpdateOn();
    rview->Update();
    viewer->redraw();
  }
}

void Fl_RViewUI::cb_applyTransformation(Fl_Button* o, void* v)
{
  rview->SetSourceTransformApply(o->value());

  // Update
  rview->SourceUpdateOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_invertTransformation(Fl_Button* o, void* v)
{
  rview->SetSourceTransformInvert(o->value());

  // Update
  rview->SourceUpdateOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_moveupTransformation(Fl_Button* o, void* v)
{
  int n;

  if (rviewUI->transformationBrowser->value() != 0) {
    n = rviewUI->transformationBrowser->value();
    if (n > 1) {
      char *buffer = strdup(rviewUI->transformationBrowser->text(n));
      rviewUI->transformationBrowser->remove(n);
      rviewUI->transformationBrowser->insert(n-1, buffer);
      rviewUI->transformationBrowser->deselect();
      rviewUI->transformationBrowser->select(n-1);
      rviewUI->info_trans_filename->value(buffer);
      free(buffer);
    }
  }
  o->clear();

  // Update
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_movedownTransformation(Fl_Button* o, void* v)
{
  int n;

  if (rviewUI->transformationBrowser->value() != 0) {
    n = rviewUI->transformationBrowser->value();
    if (n < rviewUI->transformationBrowser->size()) {
      char *buffer = strdup(rviewUI->transformationBrowser->text(n));
      rviewUI->transformationBrowser->remove(n);
      rviewUI->transformationBrowser->insert(n+1, buffer);
      rviewUI->transformationBrowser->deselect();
      rviewUI->transformationBrowser->select(n+1);
      rviewUI->info_trans_filename->value(buffer);
      free(buffer);
    }
  }
  o->clear();

  // Update
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_movieTransformation(Fl_Button* o, void* v)
{
  int i;
  char buffer2[256];

  if (fl_choice("Do you want to save the movie to disk?", NULL, "No", "Yes") == 2) {

    for (i = 1; i <= rviewUI->transformationBrowser->size(); i++) {

      // Read transformation
      char *buffer = strdup(rviewUI->transformationBrowser->text(i));
      rview->ReadTransformation(buffer);
      rviewUI->info_trans_filename->value(buffer);
      free(buffer);

      // Update
      rview->Update();
      sprintf(buffer2, "movie_%.5d.png", i);
      rview->DrawOffscreen(buffer2);
    }
  } else {

    for (i = 1; i <= rviewUI->transformationBrowser->size(); i++) {

      // Read transformation
      char *buffer = strdup(rviewUI->transformationBrowser->text(i));
      rview->ReadTransformation(buffer);
      rviewUI->info_trans_filename->value(buffer);
      free(buffer);

      // Update
      rview->Update();
      viewer->redraw();

      // Force drawing
      Fl::wait(0);
    }
  }

  if (rviewUI->transformationBrowser->size() > 0) {
    // Read transformation
    char *buffer = strdup(rviewUI->transformationBrowser->text(rviewUI->transformationBrowser->value()));
    rview->ReadTransformation(buffer);
    rviewUI->info_trans_filename->value(buffer);
    free(buffer);
  }

  // Update
  rview->Update();
  viewer->redraw();

  o->clear();
}

void Fl_RViewUI::cb_viewDeformationGridResolution(Fl_Slider* o, void* v)
{
  rview->SetDisplayDeformationGridResolution(round(o->value()));
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewDeformationGrid(Fl_Button* o, void* v)
{
  if (o->value() == 0) rview->DisplayDeformationGridOff();
  if (o->value() == 1) rview->DisplayDeformationGridOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewDeformationArrows(Fl_Button* o, void* v)
{
  if (o->value() == 0) rview->DisplayDeformationArrowsOff();
  if (o->value() == 1) rview->DisplayDeformationArrowsOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewDeformationPoints(Fl_Button* o, void* v)
{
  if (o->value() == 0) rview->DisplayDeformationPointsOff();
  if (o->value() == 1) rview->DisplayDeformationPointsOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_DeformationProperty(Fl_Menu_*, void* v)
{
  if (strcmp((char *)v, "none") == 0) {
    rview->SetDeformationProperty(None);
  }
  if (strcmp((char *)v, "displacement") == 0) {
    rview->SetDeformationProperty(Displacement);
  }
  if (strcmp((char *)v, "jacobian") == 0) {
    rview->SetDeformationProperty(Jacobian);
  }
  if (strcmp((char *)v, "jacobian_expansion") == 0) {
    rview->SetDeformationProperty(Jacobian_Expansion);
  }
  if (strcmp((char *)v, "jacobian_contraction") == 0) {
    rview->SetDeformationProperty(Jacobian_Contraction);
  }
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_deformationMinMax(Fl_Value_Slider* o, void* v)
{
  rview->GetDeformationLookupTable()->SetMinIntensity(round(rviewUI->deformationMin->value()));
  rviewUI->deformationMin->value(rview->GetDeformationLookupTable()->GetMinIntensity());
  rview->GetDeformationLookupTable()->SetMaxIntensity(round(rviewUI->deformationMax->value()));
  rviewUI->deformationMax->value(rview->GetDeformationLookupTable()->GetMaxIntensity());
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_deformationBlending(Fl_Value_Slider* o, void* v)
{
  rview->SetDeformationBlending(rviewUI->deformationBlending->value());
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::UpdateTransformationControlWindow()
{
  int j;
  list<char *> textList;
  list<char *>::iterator i;

  // Get list with transformation information
  rview->GetTransformationText(textList);

  // Clear browser
  info_trans_details->clear();

  // Add transformation information
  for (i = textList.begin(); i != textList.end(); i++) {
    info_trans_details->add(*i);
  }

  // Update transformation display
  rviewUI->transformApply->value(rview->GetSourceTransformApply());
  rviewUI->transformInvert->value(rview->GetSourceTransformInvert());
  rviewUI->viewDeformationGrid->value(rview->GetDisplayDeformationGrid());
  rviewUI->viewDeformationPoints->value(rview->GetDisplayDeformationPoints());
  rviewUI->viewDeformationArrows->value(rview->GetDisplayDeformationArrows());
  rviewUI->viewDeformationGridResolution->value(rview->GetDisplayDeformationGridResolution());

  // Get transformation
  irtkTransformation *transform = rview->GetTransformation();

  // Reset values to identity
  for (j = 0; j < transform->NumberOfDOFs(); j++) {
    rviewUI->transformationValuator[j]->value(transform->Get(j));
  }

  deformationMin->minimum(rview->GetDeformationLookupTable()->minData);
  deformationMin->maximum(rview->GetDeformationLookupTable()->maxData);
  deformationMin->value(rview->GetDeformationLookupTable()->minDisplay);
  deformationMax->minimum(rview->GetDeformationLookupTable()->minData);
  deformationMax->maximum(rview->GetDeformationLookupTable()->maxData);
  deformationMax->value(rview->GetDeformationLookupTable()->maxDisplay);
  deformationBlending->value(rview->GetDeformationBlending());
}

void Fl_RViewUI::ShowTransformationControlWindow()
{}

void Fl_RViewUI::InitializeTransformationControlWindow()
{
  {
    // Create transformation controls
    Fl_Group* o = new Fl_Group(0, 30, 400, 410, "Transformation controls");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      transformationBrowser = new Fl_Simple_Browser(10, 40, 370, 130, "Transformations");
      transformationBrowser->callback((Fl_Callback*)cb_browseTransformation);
      {
        Fl_Button  *o = new Fl_Button(15,  190,  100,  30, "Move up");
        o->callback((Fl_Callback*)cb_moveupTransformation);
      }
      {
        Fl_Button  *o = new Fl_Button(145, 190,  100,  30, "Edit");
        o->callback((Fl_Callback*)cb_editTransformation);
      }
      {
        Fl_Button  *o = new Fl_Button(275, 190,  100,  30, "Reset");
        o->callback((Fl_Callback*)cb_editTransformationReset);
      }
      {
        Fl_Button  *o = new Fl_Button(15,  230,  100,  30, "Move down");
        o->callback((Fl_Callback*)cb_movedownTransformation);
      }
      {
        Fl_Button  *o = new Fl_Button(145, 230,  100,  30, "Delete");
        o->callback((Fl_Callback*)cb_deleteTransformation);
      }
      {
        Fl_Button  *o = new Fl_Button(275, 230,  100,  30, "Movie");
        o->callback((Fl_Callback*)cb_movieTransformation);
      }

      transformApply = new Fl_Check_Button(20, 270, 170, 20, "  Apply transformation");
      transformApply->callback((Fl_Callback*)cb_applyTransformation);

      transformInvert = new Fl_Check_Button(200, 270, 170, 20, "  Invert transformation");
      transformInvert->callback((Fl_Callback*)cb_invertTransformation);

      info_trans_filename = new Fl_Output(150, 300, 230, 20, "Transform filename = ");
      info_trans_filename->box(FL_FLAT_BOX);
      info_trans_details = new Fl_Hold_Browser(150, 330, 230, 100, "Transform details    = ");
      info_trans_details->align(FL_ALIGN_LEFT);
    }
    o->end(); // End of transformation controls
  }
  {
    // Create deformation display controls
    Fl_Group* o = new Fl_Group(0, 440, 400, 220, "Deformation Display");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Check_Button *o  = viewDeformationArrows = new Fl_Check_Button(20, 455, 60, 20, "Vectors");
      o->callback((Fl_Callback*)cb_viewDeformationArrows);
    }
    {
      Fl_Check_Button *o  = viewDeformationGrid = new Fl_Check_Button(90, 455, 60, 20, "Grid");
      o->callback((Fl_Callback*)cb_viewDeformationGrid);
    }
    {
      Fl_Check_Button *o  = viewDeformationPoints = new Fl_Check_Button(140, 455, 120, 20, "Control points");
      o->callback((Fl_Callback*)cb_viewDeformationPoints);
    }
    {
      Fl_Choice* o = deformationProperty = new Fl_Choice(260, 450, 120, 30);
      o->align(FL_ALIGN_BOTTOM);
      o->menu(menu_deformationProperty);
    }
    {
      Fl_Value_Slider* o = deformationMin = new Fl_Value_Slider(10, 500, 380, 20, "Min. value");
      o->step(1);
      o->minimum(rview->GetDeformationLookupTable()->minData);
      o->maximum(rview->GetDeformationLookupTable()->maxData);
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->align(FL_ALIGN_TOP_LEFT);
      o->callback((Fl_Callback*)cb_deformationMinMax);
    }
    {
      Fl_Value_Slider* o = deformationMax = new Fl_Value_Slider(10, 540, 380, 20, "Max. value");
      o->step(1);
      o->minimum(rview->GetDeformationLookupTable()->minData);
      o->maximum(rview->GetDeformationLookupTable()->maxData);
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->align(FL_ALIGN_TOP_LEFT);
      o->callback((Fl_Callback*)cb_deformationMinMax);
    }
    {
      Fl_Value_Slider* o = deformationBlending = new Fl_Value_Slider(10, 580, 380, 20, "Blending");
      o->step(0.01);
      o->minimum(0);
      o->maximum(1);
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->align(FL_ALIGN_TOP_LEFT);
      o->callback((Fl_Callback*)cb_deformationBlending);
    }
    {
      Fl_Value_Slider* o = viewDeformationGridResolution = new Fl_Value_Slider(10, 620, 380, 20, "Deformation field resolution");
      o->step(1);
      o->minimum(0);
      o->maximum(40);
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->align(FL_ALIGN_TOP_LEFT);
      o->callback((Fl_Callback*)cb_viewDeformationGridResolution);
    }
    o->end(); // End of deformation display controls
  }
}

void Fl_RViewUI::InitializeTransformationEditor()
{
  // Global transformation Editor
  {
    editTransformation = new Fl_Window(750, 230, "editGlobalTransformation");
    editTransformation->set_modal();
    {
      Fl_Window *o = new Fl_Window(0, 0, 750, 180);
      o->box(FL_ENGRAVED_BOX);
      {
        // Boxes
        Fl_Box **o = new Fl_Box *[5];
        o[0] = new Fl_Box( 10,  40, 100, 20, "Translation (mm)");
        o[1] = new Fl_Box( 10,  70, 100, 20, "Rotation (deg)");
        o[2] = new Fl_Box( 10, 100, 100, 20, "Scaling (%)");
        o[3] = new Fl_Box( 10, 130, 100, 20, "Skews (deg)");
        o[4] = new Fl_Box(320,  10, 100, 20, "Transformation parameters");
      }
      {
        {
          // Sliders
          int i;
          Fl_Valuator** o = transformationValuator = new Fl_Valuator *[12];

          // translation
          o[0]  = new Fl_Value_Slider(140, 40, 180, 20, "tx");
          o[1]  = new Fl_Value_Slider(350, 40, 180, 20, "ty");
          o[2]  = new Fl_Value_Slider(560, 40, 180, 20, "tz");
          for (i=0; i<3; i++) {
            o[i]->range(_GLOBAL_TRANSLATION_MIN, _GLOBAL_TRANSLATION_MAX);
            o[i]->step(0.5);
          }

          // rotation
          o[3]  = new Fl_Value_Slider(140, 70, 180, 20, "rx");
          o[4]  = new Fl_Value_Slider(350, 70, 180, 20, "ry");
          o[5]  = new Fl_Value_Slider(560, 70, 180, 20, "rz");
          for (i=3; i<6; i++) {
            o[i]->range(_GLOBAL_ROTATION_MIN, _GLOBAL_ROTATION_MAX);
            o[i]->step(1);
          }

          // scaling
          o[6]  = new Fl_Value_Slider(140, 100, 180, 20, "sx");
          o[7]  = new Fl_Value_Slider(350, 100, 180, 20, "sy");
          o[8]  = new Fl_Value_Slider(560, 100, 180, 20, "sz");
          for (i=6; i<9; i++) {
            o[i]->range(_GLOBAL_SCALING_MIN, _GLOBAL_SCALING_MAX);
            o[i]->step(5);
          }

          // skews
          o[9]  = new Fl_Value_Slider(140, 130, 180, 20, "xy");
          o[10] = new Fl_Value_Slider(350, 130, 180, 20, "yz");
          o[11] = new Fl_Value_Slider(560, 130, 180, 20, "xz");
          for (i=9; i<12; i++) {
            o[i]->range(_GLOBAL_ROTATION_MIN, _GLOBAL_ROTATION_MAX);
            o[i]->step(1);
          }

          // Settings and callback
          irtkTransformation *transform = rview->GetTransformation();
          for (i=0; i<transform->NumberOfDOFs(); i++) {
            o[i]->callback((Fl_Callback*)cb_editTransformationUpdate,
                           (void *)i);
            o[i]->type(FL_HOR_NICE_SLIDER);
            o[i]->box(FL_EMBOSSED_BOX);
            o[i]->value(transform->Get(i));
            o[i]->align(FL_ALIGN_LEFT);
          }
        }
      }
      o->end();
    }
    {
      Fl_Window *o = new Fl_Window(0, 180, 750, 50);
      o->box(FL_ENGRAVED_BOX);
      {
        // Buttons
        Fl_Button **o = new  Fl_Button *[2];
        o[0] = new Fl_Return_Button( 10, 10, 80, 30, "OK");
        o[1] = new Fl_Return_Button(110, 10, 80, 30, "Reset");
        o[0]->callback((Fl_Callback *)cb_editTransformationApply);
        o[1]->callback((Fl_Callback *)cb_editTransformationReset);
      }
      o->end();
    }
    editTransformation->end(); // End of editGlobalTransformation
  }
}

#endif
