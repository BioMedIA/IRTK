/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

Copyright (c) IXICO LIMITED
All rights reserved.
See COPYRIGHT for details

=========================================================================*/

#include <irtkImage.h>
#include <irtkTransformation.h>
#include <irtkRegistration.h>

#include <sys/types.h>
#include <sys/stat.h>

// Local includes
#include <Fl_RViewUI.h>

#ifdef HAS_TRANSFORMATION_PANEL

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

  if (dynamic_cast<irtkFreeFormTransformation *>(transform) != NULL) return;

  for (i = 0; i<transform->NumberOfDOFs(); i++) {
    rviewUI->transformationValuator[i]->value(transform->Get(i));
  }

}

void Fl_RViewUI::cb_loadTransformation(Fl_Button *, void *)
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

  if (dynamic_cast<irtkFreeFormTransformation *>(transform) != NULL) return;

  for (i = 0; i<transform->NumberOfDOFs(); i++) {
    rviewUI->transformationValuator[i]->value(transform->Get(i));
  }
}

void Fl_RViewUI::cb_saveTransformation(Fl_Button *, void *)
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

void Fl_RViewUI::cb_browseTransformation(Fl_Browser* o, void*)
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

void Fl_RViewUI::cb_deleteTransformation(Fl_Button* o, void*)
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

void Fl_RViewUI::cb_editTransformation(Fl_Button* o, void*)
{
  rviewUI->editTransformation->show();
  o->clear();

  // Update
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}


void Fl_RViewUI::cb_editTransformationApply(Fl_Button*, void*)
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

void Fl_RViewUI::cb_editTransformationReset(Fl_Button*, void*)
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

void Fl_RViewUI::cb_applyTransformation(Fl_Button* o, void*)
{
  rview->SetSourceTransformApply(o->value());

  // Update
  rview->SourceUpdateOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_invertTransformation(Fl_Button* o, void*)
{
  rview->SetSourceTransformInvert(o->value());

  // Update
  rview->SourceUpdateOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_moveupTransformation(Fl_Button* o, void*)
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

void Fl_RViewUI::cb_movedownTransformation(Fl_Button* o, void*)
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

void Fl_RViewUI::cb_movieTransformation(Fl_Button* o, void*)
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

void Fl_RViewUI::cb_viewDeformationGridResolution(Fl_Slider* o, void*)
{
  rview->SetDisplayDeformationGridResolution(round(o->value()));
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewDeformationGrid(Fl_Button* o, void*)
{
  if (o->value() == 0) rview->DisplayDeformationGridOff();
  if (o->value() == 1) rview->DisplayDeformationGridOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewDeformationArrows(Fl_Button* o, void*)
{
  if (o->value() == 0) rview->DisplayDeformationArrowsOff();
  if (o->value() == 1) rview->DisplayDeformationArrowsOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewDeformationPoints(Fl_Button* o, void*)
{
  if (o->value() == 0) rview->DisplayDeformationPointsOff();
  if (o->value() == 1) rview->DisplayDeformationPointsOn();
  rview->Update();
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

  // Reset values to identity if transformation is has affine component
  if (dynamic_cast<irtkFreeFormTransformation *>(transform) != NULL) return;

  for (j = 0; j < transform->NumberOfDOFs(); j++) {
    rviewUI->transformationValuator[j]->value(transform->Get(j));
  }
}

void Fl_RViewUI::ShowTransformationControlWindow()
{}

void Fl_RViewUI::InitializeTransformationControlWindow()
{
  {
    // Create transformation controls
    Fl_Group* o = new Fl_Group(0, 30, 400, 460, "Transformation controls");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      transformationBrowser = new Fl_Simple_Browser(10, 40, 370, 180, "Transformations");
      transformationBrowser->callback((Fl_Callback*)cb_browseTransformation);
      {
        Fl_Button  *o = new Fl_Button(15,  240,  100,  30, "Move up");
        o->callback((Fl_Callback*)cb_moveupTransformation);
      }
      {
        Fl_Button  *o = new Fl_Button(145, 240,  100,  30, "Edit");
        o->callback((Fl_Callback*)cb_editTransformation);
      }
      {
        Fl_Button  *o = new Fl_Button(275, 240,  100,  30, "Reset");
        o->callback((Fl_Callback*)cb_editTransformationReset);
      }
      {
        Fl_Button  *o = new Fl_Button(15,  280,  100,  30, "Move down");
        o->callback((Fl_Callback*)cb_movedownTransformation);
      }
      {
        Fl_Button  *o = new Fl_Button(145, 280,  100,  30, "Delete");
        o->callback((Fl_Callback*)cb_deleteTransformation);
      }
      {
        Fl_Button  *o = new Fl_Button(275, 280,  100,  30, "Movie");
        o->callback((Fl_Callback*)cb_movieTransformation);
      }

      transformApply = new Fl_Check_Button(20, 320, 170, 20, "  Apply transformation");
      transformApply->callback((Fl_Callback*)cb_applyTransformation);

      transformInvert = new Fl_Check_Button(200, 320, 170, 20, "  Invert transformation");
      transformInvert->callback((Fl_Callback*)cb_invertTransformation);

      info_trans_filename = new Fl_Output(150, 350, 230, 20, "Transform filename = ");
      info_trans_filename->box(FL_FLAT_BOX);
      info_trans_details = new Fl_Hold_Browser(150, 380, 230, 100, "Transform details    = ");
      info_trans_details->align(FL_ALIGN_LEFT);
    }
    o->end(); // End of transformation controls
  }
  {
    // Create deformation display controls
    Fl_Group* o = new Fl_Group(0, 490, 400, 220, "Deformation Display");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Check_Button *o  = viewDeformationArrows = new Fl_Check_Button(20, 510, 60, 20, "Deformation vectors");
      o->callback((Fl_Callback*)cb_viewDeformationArrows);
    }
    {
      Fl_Check_Button *o  = viewDeformationGrid = new Fl_Check_Button(210, 510, 60, 20, "Deformation grid");
      o->callback((Fl_Callback*)cb_viewDeformationGrid);
    }
    {
      Fl_Check_Button *o  = viewDeformationPoints = new Fl_Check_Button(20, 550, 120, 20, "Control points");
      o->callback((Fl_Callback*)cb_viewDeformationPoints);
    }
    {
      Fl_Value_Slider* o = viewDeformationGridResolution = new Fl_Value_Slider(10, 610, 380, 20, "Deformation field resolution");
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
                           reinterpret_cast<void *>(i));
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
