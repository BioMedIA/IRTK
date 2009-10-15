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

// Local includes
#include <Fl_RViewUI.h>

extern Fl_RViewUI  *rviewUI;
extern Fl_RView    *viewer;
extern irtkRView    *rview;

#include <bitmaps/fileopen.xpm>
#include <bitmaps/filesave.xpm>
#include <bitmaps/camera.xpm>
#include <bitmaps/axial.xpm>
#include <bitmaps/coronal.xpm>
#include <bitmaps/sagital.xpm>
#include <bitmaps/multiview.xpm>
//#include <bitmaps/messagebox_info.xpm>
//#include <bitmaps/messagebox_noinfo.xpm>
#include <bitmaps/player_play.xpm>
#include <bitmaps/player_pause.xpm>
#include <bitmaps/moviecapture.xpm>

void quit(Fl_Widget*, void*)
{
  exit(0);
}

void Fl_RViewUI::cb_saveScreen(Fl_Button *, void *v)
{
  char *filename = fl_file_chooser("Save screen", "*.png", "");
  if (filename != NULL) {
    rview->DrawOffscreen(filename);
  }
}

void Fl_RViewUI::cb_ViewMode(Fl_Button *o, void *v)
{
  rview->Configure((irtkRViewConfig *)v);
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_resliceX(Fl_Roller* o, void* v)
{
  double x, y, z;

  // Calculate new reslice position
  rview->GetOrigin(x, y, z);
  rview->GetTarget()->WorldToImage(x, y, z);
  x = round(o->value());
  if (x < 0) x = 0;
  if (x >= rview->GetTarget()->GetX()) x = rview->GetTarget()->GetX()-1;
  rview->GetTarget()->ImageToWorld(x, y, z);
  rview->SetOrigin(x, y, z);

  // Update
  rview->Update();
  viewer->redraw();
  rviewUI->update();
}

void Fl_RViewUI::cb_resliceY(Fl_Roller* o, void* v)
{
  double x, y, z;

  // Calculate new reslice position
  rview->GetOrigin(x, y, z);
  rview->GetTarget()->WorldToImage(x, y, z);
  y = round(o->value());
  if (y < 0) y = 0;
  if (y >= rview->GetTarget()->GetY()) y = rview->GetTarget()->GetY()-1;
  rview->GetTarget()->ImageToWorld(x, y, z);
  rview->SetOrigin(x, y, z);

  // Update
  rview->Update();
  viewer->redraw();
  rviewUI->update();
}

void Fl_RViewUI::cb_resliceZ(Fl_Roller* o, void* v)
{
  double x, y, z;

  // Calculate new reslice position
  rview->GetOrigin(x, y, z);
  rview->GetTarget()->WorldToImage(x, y, z);
  z = round(o->value());
  if (z < 0) z = 0;
  if (z >= rview->GetTarget()->GetZ()) z = rview->GetTarget()->GetZ()-1;
  rview->GetTarget()->ImageToWorld(x, y, z);
  rview->SetOrigin(x, y, z);

  // Update
  rview->Update();
  viewer->redraw();
  rviewUI->update();
}

void Fl_RViewUI::cb_zoom(Fl_Value_Slider* o, void* v)
{
  rview->SetResolution(o->value());
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_snapGrid(Fl_Check_Button* o, void* v)
{
  if (o->value() == 1) {
    rview->SnapToGridOn();
  } else {
    rview->SnapToGridOff();
  }
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewCursor(Fl_Check_Button* o, void* v)
{
  if (o->value() == 1) {
    rview->DisplayCursorOn();
    rview->DisplayAxisLabelsOn();
  } else {
    rview->DisplayCursorOff();
    rview->DisplayAxisLabelsOff();
  }
  rview->Update();
  viewer->redraw();
}

void Fl_RViewUI::cb_viewCursorMode(Fl_Button* o, void* v)
{
  if (strcmp((char *)v, "CrossHair") == 0) {
    rview->SetCursorMode(CrossHair);
  } else if (strcmp((char *)v, "CursorX") == 0) {
    rview->SetCursorMode(CursorX);
  } else if (strcmp((char *)v, "CursorV") == 0) {
    rview->SetCursorMode(CursorV);
  } else if (strcmp((char *)v, "CursorBar") == 0) {
    rview->SetCursorMode(CursorBar);
  }
  rview->Update();
  viewer->redraw();
}

Fl_RViewUI::Fl_RViewUI()
{
  // Create main window
  Fl_Window* o = mainWindow = new Fl_Window(1060, 780, "rview");
  o->labeltype(FL_NORMAL_LABEL);
  o->box(FL_ENGRAVED_BOX);
  o->user_data((void*)(this));
  Fl_Group::current()->resizable(o);

  // Initialize the main window
  this->InitializeMainWindow();

  {
    // Create control windows
    Fl_Window* o = new Fl_Window(660, 60, 400, 660);
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    {
      // Create control tabs
      Fl_Tabs* o = tab_menu = new Fl_Tabs(0, 0, 400, 660);
      o->selection_color(15);
#ifdef HAS_VISUALISATION_PANEL
      {
        // Create image controls
        Fl_Group* o = new Fl_Group(0, 30, 400, 630, "Visualisation");
        o->user_data((void*)(this));
        // Initialize the control window
        this->InitializeImageControlWindow();
        o->end();
      }
#endif
#ifdef HAS_TRANSFORMATION_PANEL
      {
        // Create transformation controls
        Fl_Group* o = new Fl_Group(0, 30, 390, 630, "Transform/Fusion");
        o->user_data((void*)(this));
        // Initialize the control window
        this->InitializeTransformationControlWindow();
        o->end();
      }
#endif
#ifdef HAS_REGISTRATION_PANEL
      {
        // Create registration controls
        Fl_Group* o = new Fl_Group(0, 30, 390, 630, "Registration");
        o->user_data((void*)(this));
        // Initialize the control window
        this->InitializeRegistrationControlWindow();
        o->end();
      }
#endif
#ifdef HAS_LANDMARK_PANEL
      {
        // Create vtk landmarks and object controls
        Fl_Group* o = new Fl_Group(0, 30, 390, 630, "Landmarks");
        o->user_data((void*)(this));
        // Initialize the control window
        this->InitializeObjectControlWindow();
        o->end();
      }
#endif
#ifdef HAS_SEGMENTATION_PANEL
      {
        // Create Segmentation controls
        Fl_Group* o = new Fl_Group(0, 30, 390, 630, "Segmentation");
        o->user_data((void*)(this));
        // Initialize the control window
        this->InitializeSegmentationControlWindow();
        o->end();
      }
#endif
      o->end();
    }
    o->end();
  }
  {
    // Create menu button windows
    Fl_Window* o = new Fl_Window(0, 0, 1060, 60);
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    {
      Fl_Simple_Group *o = new Fl_Simple_Group(5, 5, 82, 50, "Target");
      {
        Fl_Button* o = new Fl_Button(6+5, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_loadTarget);
        o->tooltip("Load target image");
        Fl_Pixmap *p = new Fl_Pixmap(fileopen_xpm);
        p->label(o);
      }
      {
        Fl_Button* o = new Fl_Button(44+5, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_saveTarget);
        o->tooltip("Save target image");
        Fl_Pixmap *p = new Fl_Pixmap(filesave_xpm);
        p->label(o);
      }
      o->end();
    }
    {
      Fl_Simple_Group *o = new Fl_Simple_Group(92, 5, 82, 50, "Source");
      {
        Fl_Button* o = new Fl_Button(6+92, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_loadSource);
        o->tooltip("Load source image");
        Fl_Pixmap *p = new Fl_Pixmap(fileopen_xpm);
        p->label(o);
      }
      {
        Fl_Button* o = new Fl_Button(44+92, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_saveSource);
        o->tooltip("Save source image");
        Fl_Pixmap *p = new Fl_Pixmap(filesave_xpm);
        p->label(o);
      }
      o->end();
    }
    {
      Fl_Simple_Group *o = new Fl_Simple_Group(179, 5, 82, 50, "Transform");
      {
        Fl_Button* o = new Fl_Button(6+179, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_loadTransformation);
        o->tooltip("Load transformation");
        Fl_Pixmap *p = new Fl_Pixmap(fileopen_xpm);
        p->label(o);
      }
      {
        Fl_Button* o = new Fl_Button(44+179, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_saveTransformation);
        o->tooltip("Save transformation");
        Fl_Pixmap *p = new Fl_Pixmap(filesave_xpm);
        p->label(o);
      }
      o->end();
    }
    {
      Fl_Simple_Group *o = new Fl_Simple_Group(266, 5, 158, 50, "Landmarks");
      {
        Fl_Button* o = new Fl_Button(6+266, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_loadTargetLandmarks);
        o->tooltip("Load target landmarks");
        Fl_Pixmap *p = new Fl_Pixmap(fileopen_xpm);
        p->label(o);
      }
      {
        Fl_Button* o = new Fl_Button(44+266, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_saveTargetLandmarks);
        o->tooltip("Save target landmarks");
        Fl_Pixmap *p = new Fl_Pixmap(filesave_xpm);
        p->label(o);
      }
      {
        Fl_Button* o = new Fl_Button(82+266, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_loadSourceLandmarks);
        o->tooltip("Load source landmarks");
        Fl_Pixmap *p = new Fl_Pixmap(fileopen_xpm);
        p->label(o);
      }
      {
        Fl_Button* o = new Fl_Button(120+266, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_saveSourceLandmarks);
        o->tooltip("Save source landmarks");
        Fl_Pixmap *p = new Fl_Pixmap(filesave_xpm);
        p->label(o);
      }
      o->end();
    }
    {
      Fl_Simple_Group *o = new Fl_Simple_Group(429, 5, 158, 50, "Segmentation");
      {
        Fl_Button* o = new Fl_Button(6+429, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_loadSegmentation);
        o->tooltip("Load segmentation");
        Fl_Pixmap *p = new Fl_Pixmap(fileopen_xpm);
        p->label(o);
      }
      {
        Fl_Button* o = new Fl_Button(44+429, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_saveSegmentation);
        o->tooltip("Save segmentation");
        Fl_Pixmap *p = new Fl_Pixmap(filesave_xpm);
        p->label(o);
      }
      {
        Fl_Button* o = new Fl_Button(82+429, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_loadSegmentTableConfig);
        o->tooltip("Load lookup table");
        Fl_Pixmap *p = new Fl_Pixmap(fileopen_xpm);
        p->label(o);
      }
      {
        Fl_Button* o = new Fl_Button(120+429, 23, 32, 32);
        o->box(FL_NO_BOX);
        o->callback((Fl_Callback*)cb_saveSegmentTableConfig);
        o->tooltip("Save lookup table");
        Fl_Pixmap *p = new Fl_Pixmap(filesave_xpm);
        p->label(o);
      }
      o->end();
    }
    {
      Fl_Button* o = new Fl_Button(595, 5, 50, 50);
      o->box(FL_NO_BOX);
      o->callback((Fl_Callback*)cb_ViewMode, View_XY);
      o->tooltip("Axial view");
      Fl_Pixmap *p = new Fl_Pixmap(axial_xpm);
      p->label(o);
    }
    {
      Fl_Button* o = new Fl_Button(655, 5, 50, 50);
      o->box(FL_NO_BOX);
      o->callback((Fl_Callback*)cb_ViewMode, View_XZ);
      o->tooltip("Coronal view");
      Fl_Pixmap *p = new Fl_Pixmap(coronal_xpm);
      p->label(o);
    }
    {
      Fl_Button* o = new Fl_Button(715, 5, 50, 50);
      o->box(FL_NO_BOX);
      o->callback((Fl_Callback*)cb_ViewMode, View_YZ);
      o->tooltip("Sagital view");
      Fl_Pixmap *p = new Fl_Pixmap(sagital_xpm);
      p->label(o);
    }
    {
      Fl_Button* o = new Fl_Button(775, 5, 50, 50);
      o->box(FL_NO_BOX);
      o->callback((Fl_Callback*)cb_ViewMode, View_XY_XZ_YZ);
      o->tooltip("Orthogonal views");
      Fl_Pixmap *p = new Fl_Pixmap(multiview_xpm);
      p->label(o);
    }
    {
      Fl_Button* o = new Fl_Button(830, 5, 50, 50);
      o->box(FL_NO_BOX);
      o->callback((Fl_Callback*)cb_saveScreen);
      o->tooltip("Save screenshot");
      Fl_Pixmap *p = new Fl_Pixmap(camera_xpm);
      p->label(o);
    }
    {
      Fl_Button* o = new Fl_Button(880, 5, 50, 50);
      o->box(FL_NO_BOX);
      o->callback((Fl_Callback*)cb_movieStart);
      o->tooltip("Animate time frames: On");
      Fl_Pixmap *p = new Fl_Pixmap(player_play_xpm);
      p->label(o);
    }
    {
      Fl_Button* o = new Fl_Button(930, 5, 50, 50);
      o->box(FL_NO_BOX);
      o->callback((Fl_Callback*)cb_movieStop);
      o->tooltip("Animate time frames: Off");
      Fl_Pixmap *p = new Fl_Pixmap(player_pause_xpm);
      p->label(o);
    }
    {
      Fl_Button* o = new Fl_Button(980, 5, 50, 50);
      o->box(FL_NO_BOX);
      o->callback((Fl_Callback*)cb_savePlayback);
      o->tooltip("Save screenshot");
      Fl_Pixmap *p = new Fl_Pixmap(moviecapture_xpm);
      p->label(o);
    }
    o->end();
  }

  {
    // Create cursor control window
    Fl_Window* o = new Fl_Window(0, 720, 660, 60);
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    // Initialize the control window
    this->InitializeCursorControlWindow();
    o->end();
  }
  {
    // Create copyright window
    Fl_Window* o = new Fl_Window(660, 720, 400, 60);
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    Fl_Group::current()->resizable(o);
    // Initialize the control window
    {
      Fl_Output* o = new Fl_Output(50, 15, 320, 30, "Copyright");
      o->value("(C) Department of Computing, Imperial College");
      o->color(FL_GRAY);
      o->box(FL_NO_BOX);
      o->align(FL_ALIGN_CENTER);
    }
    o->end();
  }
  o->end();

#ifdef HAS_TRANSFORMATION_PANEL
  // Initialize, but don't show popup windows
  this->InitializeTransformationEditor();
#endif

  // Set up correct correct value for zoom widget
  zoom->value(viewer->v->GetResolution());
}

void Fl_RViewUI::show()
{
#ifdef HAS_VISUALISATION_PANEL
  rviewUI->ShowImageControlWindow();
#endif

#ifdef HAS_TRANSFORMATION_PANEL
  rviewUI->ShowTransformationControlWindow();
#endif

#ifdef HAS_LANDMARK_PANEL
  rviewUI->ShowObjectControlWindow();
#endif

#ifdef HAS_SEGMENTATION_PANEL
  rviewUI->ShowSegmentationControlWindow();
#endif

  Fl_Window *settingsWindow = new Fl_Window(200, 130, "Settings");
  settingsWindow->set_modal();

  /*
   * Cursor (display, shape, colour), Labels (display, colour),
   */
  {
    // Create buttons
    Fl_Button **button = new  Fl_Button *[2];
    button[0] = new Fl_Return_Button( 10, 90,  80, 30, "OK");
    button[1] = new Fl_Return_Button(110, 90,  80, 30, "Cancel");
  }
  settingsWindow->end();

  // Update viewer
  this->update();

  // Show the main window
  mainWindow->show();
}

void Fl_RViewUI::update()
{
  double x, y, z;
  char buffer1[256], buffer2[256], buffer3[256], buffer4[256], buffer5[256];

  // Update info
  viewer->v->GetInfoText(buffer1, buffer2, buffer3, buffer4, buffer5);
  info_voxel->value(buffer1);
  info_world->value(buffer2);
  info_target->value(buffer3);
  info_source->value(buffer4);
  info_segmentation->value(buffer5);
  info_snap_to_grid->value(rview->GetSnapToGrid());
  info_cursor->value(rview->GetDisplayCursor());

  // Calculate new reslice position
  rview->GetOrigin(x, y, z);
  rview->GetTarget()->WorldToImage(x, y, z);
  sliceX->value(round(x));
  sliceY->value(round(y));
  sliceZ->value(round(z));

  if (rview->GetTarget()->GetX() > 0) {
    sliceX->activate();
    sliceX->maximum(rview->GetTarget()->GetX());
  } else {
    sliceX->deactivate();
  }
  if (rview->GetTarget()->GetY() > 0) {
    sliceY->activate();
    sliceY->maximum(rview->GetTarget()->GetY());
  } else {
    sliceY->deactivate();
  }
  if (rview->GetTarget()->GetZ() > 0) {
    sliceZ->activate();
    sliceZ->maximum(rview->GetTarget()->GetZ());
  } else {
    sliceZ->deactivate();
  }

  // Update zoom
  rviewUI->zoom->value(rview->GetResolution());

  // Update panels
#ifdef HAS_VISUALISATION_PANEL
  rviewUI->UpdateImageControlWindow();
#endif
#ifdef HAS_TRANSFORMATION_PANEL
  rviewUI->UpdateTransformationControlWindow();
#endif
#ifdef HAS_REGISTRATION_PANEL
  rviewUI->UpdateRegistrationControlWindow();
#endif
#ifdef HAS_LANDMARK_PANEL
  rviewUI->UpdateObjectControlWindow();
#endif
#ifdef HAS_SEGMENTATION_PANEL
  rviewUI->UpdateSegmentationControlWindow();
#endif
}

void Fl_RViewUI::InitializeMainWindow()
{
  {
    Fl_RView* o = viewer = new Fl_RView(24, 84, 612, 612, "RViewer");
    rview = viewer->v;
    Fl_Group::current()->resizable(o);
    o->end();
  }
  {
    Fl_Value_Slider* o = zoom = new Fl_Value_Slider(190, 61, 280, 20);
    o->type(5);
    o->box(FL_EMBOSSED_BOX);
    o->maximum(10);
    o->step(0.1);
    o->value(1);
    o->callback((Fl_Callback*)cb_zoom);
    o->when(FL_WHEN_RELEASE);
  }
  {
    Fl_Roller *o = sliceX = new Fl_Roller(230, 700, 200, 20);
    o->step(1);
    o->type(FL_HORIZONTAL);
    o->deactivate();
    o->callback((Fl_Callback *)cb_resliceX);
  }
  {
    Fl_Roller *o = sliceY = new Fl_Roller(1, 290, 20, 200);
    o->step(1);
    o->deactivate();
    o->callback((Fl_Callback *)cb_resliceY);
  }
  {
    Fl_Roller *o = sliceZ = new Fl_Roller(640, 290, 20, 200);
    o->step(1);
    o->deactivate();
    o->callback((Fl_Callback *)cb_resliceZ);
  }
}

void Fl_RViewUI::InitializeCursorControlWindow()
{
  {
    Fl_Output* o = info_voxel = new Fl_Output(145, 5, 120, 20, "Voxel coordinates = ");
    o->box(FL_FLAT_BOX);
  }
  {
    Fl_Output* o = info_world = new Fl_Output(145, 30, 120, 20, "World coordinates = ");
    o->box(FL_FLAT_BOX);
  }
  {
    Fl_Output* o = info_target = new Fl_Output(350, 5, 60, 20, "Target = ");
    o->box(FL_FLAT_BOX);
  }
  {
    Fl_Output* o = info_source = new Fl_Output(480, 5, 60, 20, "Source = ");
    o->box(FL_FLAT_BOX);
  }
  {
    Fl_Output* o = info_segmentation = new Fl_Output(350, 30, 190, 20, "Anatomy = ");
    o->box(FL_FLAT_BOX);
  }
  {
    Fl_Check_Button *o = info_snap_to_grid = new Fl_Check_Button(550, 5, 100, 20, "Snap to grid");
    o->callback((Fl_Callback*)cb_snapGrid);
    o->tooltip("Snap to nearest voxel when reslicing: On/Off");
  }
  {
    Fl_Check_Button *o = info_cursor = new Fl_Check_Button(550, 30, 100, 20, "Show cursor");
    o->callback((Fl_Callback*)cb_viewCursor);
    o->tooltip("Show cursor and orientation in viewer: On/Off");
  }
}
