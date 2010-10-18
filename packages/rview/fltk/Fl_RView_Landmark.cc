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

#ifdef HAS_LANDMARK_PANEL

// Bitmaps
#include <bitmaps/landmarks.xbm>

extern Fl_RViewUI  *rviewUI;
extern Fl_RView    *viewer;
extern irtkRView    *rview;

void Fl_RViewUI::cb_loadTargetLandmarks(Fl_Button *, void *)
{
  int i;
  char *filename = fl_file_chooser("Load target landmarks", "*.vtk", "");
  if (filename != NULL) {

    // Read landmarks
    rview->ReadTargetLandmarks(filename);

    // Clear browser
    for (i = 1; i <= rviewUI->targetLandmarkBrowser->size(); i++) {
      delete (char *)(rviewUI->targetLandmarkBrowser->data(i));
    }
    rviewUI->targetLandmarkBrowser->clear();

    // Update browser
    for (i = 1; i <= rview->GetNumberOfTargetLandmarks(); i++) {
      char buffer[256], *label = NULL;
      irtkPoint point;

      // Get landmark
      rview->GetTargetLandmark(point, i, label);
      if (label == NULL) {
        label = new char [256];
        sprintf(label, "%s", ""); // To avoid zero-length format string warning
      }

      // Update browser entry
      sprintf(buffer, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
              i, label, point._x, point._y, point._z);
      rviewUI->targetLandmarkBrowser->add(buffer, (void *)label);
    }

    // Update
    rview->Update();
    viewer->redraw();
  }
}

void Fl_RViewUI::cb_loadSourceLandmarks(Fl_Button *, void *)
{
  int i;
  char *filename = fl_file_chooser("Load source landmarks", "*.vtk", "");
  if (filename != NULL) {

    // Read landmarks
    rview->ReadSourceLandmarks(filename);

    // Clear browser
    for (i = 1; i <= rviewUI->sourceLandmarkBrowser->size(); i++) {
      delete (char *)(rviewUI->sourceLandmarkBrowser->data(i));
    }
    rviewUI->sourceLandmarkBrowser->clear();

    // Update browser
    for (i = 1; i <= rview->GetNumberOfSourceLandmarks(); i++) {
      char buffer[256], *label = NULL;
      irtkPoint point;

      // Get landmark
      rview->GetSourceLandmark(point, i, label);
      if (label == NULL) {
        label = new char [256];
        sprintf(label, "%s", ""); // To avoid zero-length format string warning
      }

      // Update browser entry
      sprintf(buffer, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
              i, label, point._x, point._y, point._z);
      rviewUI->sourceLandmarkBrowser->add(buffer, (void *)label);
    }

    // Update
    rview->Update();
    viewer->redraw();
  }
}

void Fl_RViewUI::cb_saveTargetLandmarks(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Save target landmarks", "*.vtk", "");
  if (filename != NULL) {
    struct stat buffer;
    if (stat(filename, &buffer) == 0) {
      if (fl_choice("File already exists. Do you want to overwrite the file?", NULL, "No", "Yes") == 2) {
        rview->WriteTargetLandmarks(filename);
      }
    } else {
      rview->WriteTargetLandmarks(filename);
    }
  }
}

void Fl_RViewUI::cb_saveSourceLandmarks(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Save source landmarks", "*.vtk", "");
  if (filename != NULL) {
    struct stat buffer;
    if (stat(filename, &buffer) == 0) {
      if (fl_choice("File already exists. Do you want to overwrite the file?", NULL, "No", "Yes") == 2) {
        rview->WriteSourceLandmarks(filename);
      }
    } else {
      rview->WriteSourceLandmarks(filename);
    }
  }
}

void Fl_RViewUI::cb_addLandmark(Fl_Button *, void* v)
{
  char buffer[256], *label = NULL;
  irtkPoint point;

  // Set landmark at cursor origin
  rview->GetOrigin(point._x, point._y, point._z);

  // Add landmark
  label = new char[256];
  sprintf(label, "%s", ""); // To avoid zero-length format string warning
  if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
    /*temp remove  if (rview->GetTrackTAG()){
	  ((irtkGenericImage<irtkGreyPixel>*)(rview->GetTarget()))->GetMaxPosition(point,1);
	  ((irtkGenericImage<irtkGreyPixel>*)(rview->GetTarget()))->GetWindowCenter(point,2);
    }*/
    rview->AddTargetLandmark(point, label);
  } else {
    /*temp remove  if (rview->GetTrackTAG()){
	  ((irtkGenericImage<irtkGreyPixel>*)(rview->GetSource()))->GetMaxPosition(point,1);
	  ((irtkGenericImage<irtkGreyPixel>*)(rview->GetSource()))->GetWindowCenter(point,2);
    }*/
    rview->AddSourceLandmark(point, label);
  }

  // Update browser
  sprintf(buffer, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
          ((Fl_Browser *)v)->size()+1, label,
          point._x, point._y, point._z);
  ((Fl_Browser *)v)->add(buffer, (void *)label);

  // Update
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_deleteLandmark(Fl_Button *, void* v)
{
  int i, id;

  // Get landmark id
  id = ((Fl_Browser *)v)->value();
  if ((id <= 0) || (id > ((Fl_Browser *)v)->size())) {
    // Alert user
    fl_alert("Need to select landmark to delete\n");
    return;
  }

  // Delete landmark
  if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
    rview->DeleteTargetLandmark(id);
  } else {
    rview->DeleteSourceLandmark(id);
  }

  // Update browser
  delete (char *)(((Fl_Browser *)v)->data(id));
  ((Fl_Browser *)v)->remove(id);

  // Update all remaining browser entries
  for (i = id; i <= ((Fl_Browser *)v)->size(); i++) {
    char buffer[256], *label = NULL;
    irtkPoint point;

    // Get landmark
    if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
      rview->GetTargetLandmark(point, i, label);
    } else {
      rview->GetSourceLandmark(point, i, label);
    }
    if (label == NULL) {
      label = (char *)(((Fl_Browser *)v)->data(i));
    }

    // Update browser entry
    sprintf(buffer, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
            i, label, point._x, point._y, point._z);
    ((Fl_Browser *)v)->remove(i);
    ((Fl_Browser *)v)->insert(i, buffer, (void *)label);
  }

  // Update
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_toggleLandmark(Fl_Input *, void* v)
{
  if (((Fl_Browser *)v)->size() == 0) {
    return;
  }

  int id;
  char *label = NULL;
  irtkPoint point;

  // Go to next landmark position
  id = ((Fl_Browser *)v)->value();   // get current landmark id
  ((Fl_Browser *)v)->deselect();     // deselect all lines
  ((Fl_Browser *)v)->value((id + 1) % ((Fl_Browser *)v)->size()); // next id
  id = ((Fl_Browser *)v)->value();   // get next landmark id
  ((Fl_Browser *)v)->middleline(id); // scroll browser
  if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
    rview->GetTargetLandmark(point, id, label);
  } else {
    rview->GetSourceLandmark(point, id, label);
  }

  // Set origin to landmark
  rview->SetOrigin(point._x, point._y, point._z);

  // Update
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_insertLandmark(Fl_Button *, void* v)
{
  int id, i;
  char buffer[256], *label = NULL;
  irtkPoint point;

  // Get landmark id
  id = ((Fl_Browser *)v)->value();
  if ((id <= 0) || (id > ((Fl_Browser *)v)->size())) {
    // Alert user
    fl_alert("Need to select landmark before which to insert\n");
    return;
  }

  // Set landmark at cursor origin
  rview->GetOrigin(point._x, point._y, point._z);

  // Insert landmark
  label = new char[256];
  sprintf(label, "%s", ""); // To avoid zero-length format string warning
  if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
    if (rview->GetTrackTAG()){
	  /*temp remove  ((irtkGenericImage<irtkGreyPixel>*)(rview->GetTarget()))->GetMaxPosition(point,1);
	  ((irtkGenericImage<irtkGreyPixel>*)(rview->GetTarget()))->GetWindowCenter(point,2);
	  */
    }
    rview->InsertTargetLandmark(point, id, label);
  } else {
    if (rview->GetTrackTAG()){
	  /*temp remove ((irtkGenericImage<irtkGreyPixel>*)(rview->GetSource()))->GetMaxPosition(point,1);
	  ((irtkGenericImage<irtkGreyPixel>*)(rview->GetSource()))->GetWindowCenter(point,2);
	  */
    }
    rview->InsertSourceLandmark(point, id, label);
  }

  // Update browser
  sprintf(buffer, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
          id, label, point._x, point._y, point._z);
  ((Fl_Browser *)v)->insert(id, buffer, (void *)label);

  // Update all following browser entries
  for (i = id+1; i <= ((Fl_Browser *)v)->size(); i++) {

    // Get landmark
    char *label2 = NULL;
    if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
      rview->GetTargetLandmark(point, i, label2);
    } else {
      rview->GetSourceLandmark(point, i, label2);
    }
    if (label2 == NULL) {
      label2 = (char *)(((Fl_Browser *)v)->data(i));
    }

    // Update browser entry
    sprintf(buffer, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
            i, label2, point._x, point._y, point._z);
    ((Fl_Browser *)v)->remove(i);
    ((Fl_Browser *)v)->insert(i, buffer, label2);
  }

  // Update
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_replaceLandmark(Fl_Button *, void* v)
{
  int id;
  char buffer[256], *label = NULL;
  irtkPoint point;

  // Get landmark id
  id = ((Fl_Browser *)v)->value();
  if ((id <= 0) || (id > ((Fl_Browser *)v)->size())) {
    // Alert user
    fl_alert("Need to select landmark to replace\n");
    return;
  }

  // Set landmark at cursor origin
  rview->GetOrigin(point._x, point._y, point._z);

  // Replace landmark
  if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
    rview->PutTargetLandmark(point, id, label);
  } else {
    rview->PutSourceLandmark(point, id, label);
  }
  label = (char *)(((Fl_Browser *)v)->data(id));

  // Update browser
  sprintf(buffer, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
          id, label, point._x, point._y, point._z);
  ((Fl_Browser *)v)->remove(id);
  ((Fl_Browser *)v)->insert(id, buffer, (void *)label);

  // Update
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_editLandmark(Fl_Button *, void* v)
{
  // Alert user
  fl_alert("Edit landmark not yet implemented");

  int id;
  const char *buffer = NULL;
  char *label = NULL;
  irtkPoint point;

  // Get landmark id
  id = ((Fl_Browser *)v)->value();
  if ((id <= 0) || (id > ((Fl_Browser *)v)->size())) {
    // Alert user
    fl_alert("Need to select landmark to label\n");
    return;
  }

  // Get landmark
  if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
    rview->GetTargetLandmark(point, id, label);
  } else {
    rview->GetSourceLandmark(point, id, label);
  }
  if (label == NULL) {
    label = (char *)(((Fl_Browser *)v)->data(id));
  }

  // Read in new label
  buffer = fl_input("Label", label);

  // Set new label
  if (buffer != NULL) {
    char text[256];

    // Copy buffer into label
    delete label;
    label = new char [256];
    sprintf(label, "%s", buffer);

    // Label landmark (does nothing at the moment)
    if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
      rview->LabelTargetLandmark(id, label);
    } else {
      rview->LabelSourceLandmark(id, label);
    }

    // Update browser
    sprintf(text, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
            id, label, point._x, point._y, point._z);
    ((Fl_Browser *)v)->text(id, text);
    ((Fl_Browser *)v)->data(id, (void *)label);
  }

  // No further update needed
}

void Fl_RViewUI::cb_browseLandmark(Fl_Browser* o, void* v)
{
  // Go to selected landmark position, if mouse button was pressed twice
  if (Fl::event_clicks() != 0) {
    char *label = NULL;
    irtkPoint point;

    // Get landmark
    if ((Fl_Browser *)v == rviewUI->targetLandmarkBrowser) {
      rview->GetTargetLandmark(point, o->value(), label);
    } else {
      rview->GetSourceLandmark(point, o->value(), label);
    }

    // Set origin to landmark
    rview->SetOrigin(point._x, point._y, point._z);

    // Update
    rview->Update();
    rviewUI->update();
    viewer->redraw();

  } else if (Fl::event_is_click() == 0) {
    // Landmark is being dragged
    cerr << "Dragging landmark not implemented, ignoring...\n";
  }
}

void Fl_RViewUI::cb_viewLandmarks(Fl_Button* o, void*)
{
  if (o->value() == 0) rview->DisplayLandmarksOff();
  if (o->value() == 1) rview->DisplayLandmarksOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_fitLandmarks(Fl_Button *, void*)
{
  if (rview->GetNumberOfTargetLandmarks() !=
      rview->GetNumberOfSourceLandmarks()) {
    fl_alert("Number of target and source landmarks does not correspond!");
    return;
  } else {
    char buffer[256];
    double error;
    error = rview->FitLandmarks();
    rview->SourceUpdateOn();
    rview->Update();
    viewer->redraw();
    sprintf(buffer, "Residual fitting error: %f mm", error);
    fl_alert(buffer);
  }
}

void Fl_RViewUI::cb_viewROI(Fl_Button* o, void*)
{
  if (o->value() == 0) rview->DisplayROIOff();
  if (o->value() == 1) rview->DisplayROIOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_trackTAG(Fl_Button* o, void*)
{
  if (o->value() == 0) rview->TrackTAGOff();
  if (o->value() == 1) rview->TrackTAGOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewTagGrid(Fl_Button* o, void*)
{
  if (o->value() == 0) rview->ViewTAGOff();
  if (o->value() == 1) rview->ViewTAGOn();
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::ShowObjectControlWindow()
{
  int i;
  char buffer[256];

  // Set up correct values for object widgets
  for (i = 1; i <= rview->GetNumberOfTargetLandmarks(); i++) {
    // Get landmark
    char *label = NULL;
    irtkPoint point;
    rview->GetTargetLandmark(point, i, label);
    if (label == NULL) {
      label = new char [256];
      sprintf(label, "%s", ""); // To avoid zero-length format string warning
    }

    // Update browser entry
    sprintf(buffer, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
            i, label, point._x, point._y, point._z);
    targetLandmarkBrowser->add(buffer, (void *)label);
  }
  for (i = 1; i <= rview->GetNumberOfSourceLandmarks(); i++) {
    // Get landmark
    char *label = NULL;
    irtkPoint point;
    rview->GetSourceLandmark(point, i, label);
    if (label == NULL) {
      label = new char [256];
      sprintf(label, "%s", ""); // To avoid zero-length format string warning
    }

    // Update browser entry
    sprintf(buffer, "%d: %s (%.2f mm, %.2f mm, %.2f mm)",
            i, label, point._x, point._y, point._z);
    sourceLandmarkBrowser->add(buffer, (void *)label);
  }
}

void Fl_RViewUI::UpdateObjectControlWindow()
{
  rviewUI->viewLandmarks->value(rview->GetDisplayLandmarks());
  rviewUI->viewTagGrid->value(rview->GetViewTAG());
}

void Fl_RViewUI::InitializeObjectControlWindow()
{
  {
    // Create target landmark controls
    Fl_Group* o = new Fl_Group(0, 30, 400, 235, "Target landmarks");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Browser *o1 = targetLandmarkBrowser =
                         new Fl_Multi_Browser(10, 40, 370, 120, "Target landmarks");
      Fl_Button  *o2 = new Fl_Return_Button(15,  180,  100,  30, "Add");
      Fl_Button  *o3 = new Fl_Return_Button(145, 180,  100,  30, "Delete");
      Fl_Button  *o4 = new Fl_Return_Button(275, 180,  100,  30, "Toggle");
      Fl_Button  *o5 = new Fl_Return_Button(15,  220,  100,  30, "Insert");
      Fl_Button  *o6 = new Fl_Return_Button(145, 220,  100,  30, "Replace");
      Fl_Button  *o7 = new Fl_Return_Button(275, 220,  100,  30, "Edit");
      o1->callback((Fl_Callback*)cb_browseLandmark, (void *)o1);
      o2->callback((Fl_Callback*)cb_addLandmark,    (void *)o1);
      o3->callback((Fl_Callback*)cb_deleteLandmark, (void *)o1);
      o4->callback((Fl_Callback*)cb_toggleLandmark, (void *)o1);
      o5->callback((Fl_Callback*)cb_insertLandmark, (void *)o1);
      o6->callback((Fl_Callback*)cb_replaceLandmark,(void *)o1);
      o7->callback((Fl_Callback*)cb_editLandmark,   (void *)o1);
    }
    o->end(); // End of target landmark controls
  }
  {
    // Create source landmark controls
    Fl_Group* o = new Fl_Group(0, 265, 400, 235, "Source landmarks");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Browser *o1 = sourceLandmarkBrowser =
                         new Fl_Multi_Browser(10,  275, 370, 120, "Source landmarks");
      Fl_Button  *o2 = new Fl_Return_Button(15,  415,  100,  30, "Add");
      Fl_Button  *o3 = new Fl_Return_Button(145, 415,  100,  30, "Delete");
      Fl_Button  *o4 = new Fl_Return_Button(275, 415,  100,  30, "Toggle");
      Fl_Button  *o5 = new Fl_Return_Button(15,  455,  100,  30, "Insert");
      Fl_Button  *o6 = new Fl_Return_Button(145, 455,  100,  30, "Replace");
      Fl_Button  *o7 = new Fl_Return_Button(275, 455,  100,  30, "Edit");
      o1->callback((Fl_Callback*)cb_browseLandmark, (void *)o1);
      o2->callback((Fl_Callback*)cb_addLandmark,    (void *)o1);
      o3->callback((Fl_Callback*)cb_deleteLandmark, (void *)o1);
      o4->callback((Fl_Callback*)cb_toggleLandmark, (void *)o1);
      o5->callback((Fl_Callback*)cb_insertLandmark, (void *)o1);
      o6->callback((Fl_Callback*)cb_replaceLandmark,(void *)o1);
      o7->callback((Fl_Callback*)cb_editLandmark,   (void *)o1);
    }
    o->end(); // End of source landmark controls
  }
  {
    // Create object display and fitting controls
    Fl_Group* o = new Fl_Group(0, 500, 400, 160, "Object Display");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Button* o = viewLandmarks = new Fl_Button(39, 520, 52, 52);
      o->type(FL_TOGGLE_BUTTON);
      o->callback((Fl_Callback*)cb_viewLandmarks);
      Fl_Bitmap *p = new Fl_Bitmap(landmarks_bits, 50, 50);
      p->label(o);
    }
    {
      Fl_Button* o = new Fl_Return_Button(125, 520, 52, 52, "FIT");
      o->callback((Fl_Callback*)cb_fitLandmarks);
    }
    {
      Fl_Button* o = new Fl_Button(211, 520, 52, 52, "ROI");
      o->callback((Fl_Callback*)cb_viewROI);
      o->type(FL_TOGGLE_BUTTON);
    }
	{
      Fl_Button* o = new Fl_Button(297, 520, 52, 52, "TAG");
      o->callback((Fl_Callback*)cb_trackTAG);
      o->type(FL_TOGGLE_BUTTON);
    }
	{
      Fl_Check_Button *o  = viewTagGrid = new Fl_Check_Button(39, 580, 60, 20, "Tag grid");
      o->callback((Fl_Callback*)cb_viewTagGrid);
    }
    o->end(); // End of object display controls
  }
}

#endif
