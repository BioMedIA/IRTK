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
#include <irtkSegment.h>

#include <FL/Fl_Window.H>

// Local includes
#include <Fl_RViewUI.h>

#ifdef HAS_SEGMENTATION_PANEL

#include <bitmaps/kchart.xpm>
#include <bitmaps/color_line_closed.xpm>
#include <bitmaps/color_line_open.xpm>

extern Fl_RViewUI  *rviewUI;
extern Fl_RView    *viewer;
extern irtkRView    *rview;

char widthStrings[10][10] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};

Fl_Menu_Item Fl_RViewUI::menu_PaintBrushWidth[] = {
  {"1", 0, (Fl_Callback*)cb_SetPaintBrushWidth, widthStrings[1]},
  {"2", 0, (Fl_Callback*)cb_SetPaintBrushWidth, widthStrings[2]},
  {"3", 0, (Fl_Callback*)cb_SetPaintBrushWidth, widthStrings[3]},
  {"4", 0, (Fl_Callback*)cb_SetPaintBrushWidth, widthStrings[4]},
  {"5", 0, (Fl_Callback*)cb_SetPaintBrushWidth, widthStrings[5]},
  {"6", 0, (Fl_Callback*)cb_SetPaintBrushWidth, widthStrings[6]},
  {"7", 0, (Fl_Callback*)cb_SetPaintBrushWidth, widthStrings[7]},
  {"8", 0, (Fl_Callback*)cb_SetPaintBrushWidth, widthStrings[8]},
  {"9", 0, (Fl_Callback*)cb_SetPaintBrushWidth, widthStrings[9]},
  {0}
};

char growStrings[2][10] = {"2D", "3D"};

Fl_Menu_Item Fl_RViewUI::menu_RegionGrowing[] = {
  {"2D", 0, (Fl_Callback*)cb_SetRegionGrowingMode, growStrings[0]},
  {"3D", 0, (Fl_Callback*)cb_SetRegionGrowingMode, growStrings[1]},
  {0}
};

void Fl_RViewUI::cb_loadSegmentation(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Load segmentation image", "*.{gipl,gipl.gz,hdr,hdr.gz,nii,nii.gz,vtk}", "");
  if (filename != NULL) {
    rview->ReadSegmentation(filename);
    //rviewUI->_id = -1;
    //rviewUI->_selected = 0;
    if (rviewUI->_histogramWindow!=NULL) rviewUI->_histogramWindow->recalculate();

    // Update
    rview->SegmentationUpdateOn();
    rview->Update();
    viewer->redraw();
    rviewUI->update();
  }
}

void Fl_RViewUI::cb_saveSegmentation(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Save segmentation image", "*.{gipl,gipl.gz,hdr,hdr.gz,nii,nii.gz,vtk}", "");
  if (filename != NULL) {
    rview->WriteSegmentation(filename);
  }
}

void Fl_RViewUI::cb_resetSegmentation(Fl_Button *, void *)
{
}

void Fl_RViewUI::cb_saveSegmentTableConfig(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Save object lookup table", "*.seg", "");

  if (filename != NULL) {
    rview->GetSegmentTable()->Write(filename);
  }
}

void Fl_RViewUI::cb_loadSegmentTableConfig(Fl_Button *, void *)
{
  char *filename = fl_file_chooser("Load object lookup table", "*.seg", "");

  if (filename != NULL) {
    rview->GetSegmentTable()->Read(filename);
    rviewUI->_id = -1;
    rviewUI->_selected = 0;
    if (rviewUI->_histogramWindow!=NULL) rviewUI->_histogramWindow->recalculate();

    // Update
    rviewUI->UpdateSegmentationBrowser();
    rview->SegmentationUpdateOn();
    rview->Update();
    rviewUI->update();
    viewer->redraw();

  }
}

void Fl_RViewUI::cb_resetSegmentTableConfig(Fl_Button *, void *)
{
  rview->GetSegmentTable()->Clear();
  rviewUI->_id = -1;
  rviewUI->_selected = 0;

  // Update
  rviewUI->UpdateSegmentationBrowser();
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_browseSegmentObject(Fl_Browser *o, void *)
{
  int id = o->value();
  //cerr <<id <<endl;

  if (id == 0) {
    rviewUI->_id = -1;
    rviewUI->_selected = 0;
  } else {
    sscanf(rviewUI->segmentObjectBrowser->text(id), "%d %*s",  &rviewUI->_id);
    sprintf(rviewUI->_label, "%s", rview->GetSegmentTable()->Get(rviewUI->_id, &rviewUI->_red, &rviewUI->_green, &rviewUI->_blue, &rviewUI->_trans, &rviewUI->_vis));
    rviewUI->_selected = o->selected(id);
    o->deselect();
    if (rviewUI->_selected) o->select(id);

  }

  // Update
  rviewUI->UpdateSegmentationControlWindow();
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_addSegment(Fl_Button *, void *)
{
  int i;
  char buffer[256];

  // Compute default id
  for (i = 0; i <= rview->GetSegmentTable()->Size(); i++) {
    if (rview->GetSegmentTable()->IsValid(i) != true) {
      rviewUI->_id = i;
      break;
    }
  }

  // Put default id
  sprintf(buffer, "%d", rviewUI->_id);

  // Ask for id
  const char *tmp = fl_input("Label ID", buffer);
  if (tmp == NULL) return;
  rviewUI->AddSegment(atoi(tmp));
}

void Fl_RViewUI::AddSegment(int label)
{
  if (rview->GetSegmentTable()->IsValid(label) == true) {
    // Label already exists
    fl_alert("Label ID already exists. Can't create label ID.");
    return;
  } else {
    // Label doesn't exist, so everything is fine
    rviewUI->_id = label;
  }

  sprintf(rviewUI->_label, "Default");
  rviewUI->_red   = 255;
  rviewUI->_green = 255;
  rviewUI->_blue  = 255;
  rviewUI->_trans = 1;
  rviewUI->_vis   = true;

  rviewUI->editSegmentLabel->value(rviewUI->_label);
  rviewUI->editSegmentTransperancy->value(rviewUI->_trans);
  rviewUI->editSegmentVisibility->value(rviewUI->_vis);
  rviewUI->editSegmentColor->color(fl_rgb_color(rviewUI->_red, rviewUI->_green, rviewUI->_blue));
  rview->GetSegmentTable()->Set(rviewUI->_id, rviewUI->_label, rviewUI->_red, rviewUI->_green, rviewUI->_blue, rviewUI->_trans, rviewUI->_vis);

  // Update
  rviewUI->UpdateSegmentationBrowser();
  rviewUI->UpdateSegmentationControlWindow();
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_deleteSegment(Fl_Button *, void *)
{
  if (rviewUI->_id == -1) {
    cerr << "ID = -1, this should never happen" << endl;
  } else {
    rview->GetSegmentTable()->Clear(rviewUI->_id);
  }
  rviewUI->_id = -1;
  rviewUI->_selected = 0;

  // Update
  rviewUI->UpdateSegmentationBrowser();
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_displaySegmentLabels(Fl_Check_Button *o, void *)
{
  if (o->value() == 1) {
    rview->SegmentationLabelsOn();
  } else {
    rview->SegmentationLabelsOff();
  }
  // Update
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_displaySegmentContours(Fl_Check_Button *o, void *)
{
  if (o->value() == 1) {
    rview->SegmentationContoursOn();
  } else {
    rview->SegmentationContoursOff();
  }
  // Update
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_selectAll(Fl_Check_Button *, void *)
{
  //if (o->value() == 1){
  for (int j=0; j<100; j++)
    if ( rview->GetSegmentTable()->IsValid(j) == true) {
      rview->GetSegmentTable()->SetVisibility(j,1);
    }
  //}
  // Update
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_deselectAll(Fl_Check_Button *, void *)
{
  for (int j=0; j<100; j++) {
    if (rview->GetSegmentTable()->IsValid(j) == true) {
      rview->GetSegmentTable()->SetVisibility(j, 0);
    }
  }
  // Update
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_editSegmentLabel(Fl_Input* o, void *)
{
  sprintf(rviewUI->_label, "%s", o->value());

  if (rviewUI->_id == -1) {
    cerr << "ID = -1, this should never happen" << endl;
  } else {
    rview->GetSegmentTable()->Set(rviewUI->_id, rviewUI->_label, rviewUI->_red, rviewUI->_green, rviewUI->_blue, rviewUI->_trans, rviewUI->_vis);
  }

  // Update
  rviewUI->UpdateSegmentationBrowser();
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_editSegmentPickColor(Fl_Button *, void *)
{
  double r = rviewUI->_red / 255.0;
  double g = rviewUI->_green / 255.0;
  double b = rviewUI->_blue / 255.0;
  int color = fl_color_chooser("Color Selection", r, g, b);

  if (color > 0) {
    rviewUI->_red   = round(r * 255);
    rviewUI->_green = round(g * 255);
    rviewUI->_blue  = round(b * 255);
  }

  // Update color
  rviewUI->editSegmentColor->color(fl_rgb_color(rviewUI->_red, rviewUI->_green, rviewUI->_blue));

  if (rviewUI->_id == -1) {
    cerr << "ID = -1, this should never happen" << endl;
  } else {
    rview->GetSegmentTable()->Set(rviewUI->_id, rviewUI->_label, rviewUI->_red, rviewUI->_green, rviewUI->_blue, rviewUI->_trans, rviewUI->_vis);
  }
  // Update
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_editSegmentTransperancy(Fl_Slider* o, void *)
{
  rviewUI->_trans = o->value();

  if (rviewUI->_id == -1) {
    cerr << "ID = -1, this should never happen" << endl;
  } else {
    rview->GetSegmentTable()->Set(rviewUI->_id, rviewUI->_label, rviewUI->_red, rviewUI->_green, rviewUI->_blue, rviewUI->_trans, rviewUI->_vis);
  }
  // Update
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_editSegmentVisibility(Fl_Check_Button* o, void *)
{
  rviewUI->_vis = int(o->value());

  if (rviewUI->_id == -1) {
    cerr << "ID = -1, this should never happen" << endl;
  } else {
    rview->GetSegmentTable()->Set(rviewUI->_id, rviewUI->_label, rviewUI->_red, rviewUI->_green, rviewUI->_blue, rviewUI->_trans, rviewUI->_vis);
  }
  // Update
  rview->SegmentationUpdateOn();
  rview->Update();
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_showHistogram(Fl_Button *, void *)
{
  if (rviewUI->_histogramWindow == NULL) {
    rviewUI->_histogramWindow = new Fl_HistogramWindow(100, 100, 500, 300, "Histogram", rview);
    Fl_Group::current()->resizable(rviewUI->_histogramWindow);
    rviewUI->_histogramWindow->color(FL_GRAY);
    rviewUI->_histogramWindow->end();
  }
  rviewUI->_histogramWindow->recalculate();
  rviewUI->_histogramWindow->show();
}

void Fl_RViewUI::cb_DrawContour(Fl_Button *, void *)
{
  rview->SegmentationMode(1);
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_DrawPaintBrush(Fl_Button *, void *)
{
  rview->SegmentationMode(0);
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_SetRegionGrowingMode(Fl_Button *, void *v)
{
  if (strcmp((char *)v, "2D") == 0) {
    rview->SetRegionGrowingMode(RegionGrowing2D);
  } else {
    rview->SetRegionGrowingMode(RegionGrowing3D);
  }
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_SetPaintBrushWidth(Fl_Button *, void *v)
{
  rview->SetPaintBrushWidth(atoi((char *)v));
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_SetRegionGrowingThresholdMinimum(Fl_Value_Slider* o, void *)
{
  rview->SetRegionGrowingThresholdMinimum(round(o->value()));
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::cb_SetRegionGrowingThresholdMaximum(Fl_Value_Slider* o, void *)
{
  rview->SetRegionGrowingThresholdMaximum(round(o->value()));
  rviewUI->update();
  viewer->redraw();
}

void Fl_RViewUI::ShowSegmentationControlWindow()
{
}

void Fl_RViewUI::UpdateSegmentationBrowser()
{
  int i;

  // Clear browser
  rviewUI->segmentObjectBrowser->clear();

  // Add labels
  for (i = 0; i <= rview->GetSegmentTable()->Size(); i++) {
    if (rview->GetSegmentTable()->IsValid(i) == true) {
      char buffer[256];
      sprintf(buffer, "%d \t %s", i, rview->GetSegmentTable()->GetLabel(i));
      rviewUI->segmentObjectBrowser->add(buffer);
    }
  }
}

void Fl_RViewUI::UpdateSegmentationControlWindow()
{
  regionGrowingMin->minimum(rview->GetTargetMin());
  regionGrowingMin->maximum(rview->GetTargetMax());
  regionGrowingMax->minimum(rview->GetTargetMin());
  regionGrowingMax->maximum(rview->GetTargetMax());
  regionGrowingMin->value(rview->GetRegionGrowingThresholdMinimum());
  regionGrowingMax->value(rview->GetRegionGrowingThresholdMaximum());

  if (rviewUI->_id == -1) {
    rviewUI->editSegmentColor->color(fl_rgb_color(255, 255, 255));
    rviewUI->editSegmentColor->deactivate();
    rviewUI->editSegmentTransperancy->value(1);
    rviewUI->editSegmentTransperancy->deactivate();
    rviewUI->editSegmentVisibility->value(1);
    rviewUI->editSegmentVisibility->deactivate();
    rviewUI->editSegmentLabel->value("");
    rviewUI->editSegmentLabel->deactivate();
  } else {
    rviewUI->editSegmentColor->color(fl_rgb_color(rviewUI->_red, rviewUI->_green, rviewUI->_blue));
    rviewUI->editSegmentColor->activate();
    rviewUI->editSegmentTransperancy->value(rviewUI->_trans);
    rviewUI->editSegmentTransperancy->activate();
    rviewUI->editSegmentVisibility->value(rviewUI->_vis);
    rviewUI->editSegmentVisibility->activate();
    rviewUI->editSegmentLabel->value(rviewUI->_label);
    rviewUI->editSegmentLabel->activate();
  }
  rviewUI->editSegmentColor->redraw();

  // Update viewing mode
  rviewUI->drawContour->value(0);
  rviewUI->drawPaintBrush->value(0);
  switch (rview->GetSegmentationMode()) {
  case 0:
    rviewUI->drawPaintBrush->value(1);
    break;
  case 1:
    rviewUI->drawContour->value(1);
    break;
  }

  switch (rview->GetRegionGrowingMode()) {
  case RegionGrowing2D:
    rviewUI->regionGrowingMode->value(0);
    break;
  case RegionGrowing3D:
    rviewUI->regionGrowingMode->value(1);
    break;
  }

  switch (rview->GetPaintBrushWidth()) {
  case 1:
    rviewUI->paintBrushWidth->value(0);
    break;
  case 2:
    rviewUI->paintBrushWidth->value(1);
    break;
  case 3:
    rviewUI->paintBrushWidth->value(2);
    break;
  case 4:
    rviewUI->paintBrushWidth->value(3);
    break;
  case 5:
    rviewUI->paintBrushWidth->value(4);
    break;
  case 6:
    rviewUI->paintBrushWidth->value(5);
    break;
  case 7:
    rviewUI->paintBrushWidth->value(6);
    break;
  case 8:
    rviewUI->paintBrushWidth->value(7);
    break;
  case 9:
    rviewUI->paintBrushWidth->value(8);
    break;
  }

  if ((rviewUI->_histogramWindow != NULL) && (rviewUI->_histogramWindow->shown())) rviewUI->_histogramWindow->draw();
}

void Fl_RViewUI::InitializeSegmentationControlWindow()
{
  _id = -1;
  _selected=0;
  {
    // Create SegmentBrowserWindow controls
    Fl_Group* o = new Fl_Group(0, 30, 400, 240, "SegmentBrowserWindow");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    o->labeltype(FL_NO_LABEL);
    {
      Fl_Multi_Browser *o1 = segmentObjectBrowser =
                               new Fl_Multi_Browser(10, 40, 370, 140);
      Fl_Button  *o2 = new Fl_Button( 15, 190,  100,  30, "Add");
      Fl_Button  *o3 = new Fl_Button(145, 190,  100,  30, "Delete");
      Fl_Check_Button  *o4 = new Fl_Check_Button(275, 190,  100,  30, "Labels");
      Fl_Check_Button  *o5 = new Fl_Check_Button(275, 230,  100,  30, "Contours");
      Fl_Button  *o6 = new Fl_Button( 15, 230,  100,  30, "Deselect all");
      Fl_Button  *o7 = new Fl_Button(145, 230,  100,  30, "Select all");
      o1->callback((Fl_Callback*)cb_browseSegmentObject);
      o2->callback((Fl_Callback*)cb_addSegment);
      o3->callback((Fl_Callback*)cb_deleteSegment);
      o4->callback((Fl_Callback*)cb_displaySegmentLabels);
      o5->callback((Fl_Callback*)cb_displaySegmentContours);
      o6->callback((Fl_Callback*)cb_deselectAll);
      o7->callback((Fl_Callback*)cb_selectAll);
      o4->align(FL_ALIGN_WRAP);
      o5->align(FL_ALIGN_WRAP);
      displayLabels = o4;
    }
    o->end(); // End of SegmentBrowserWindow controls

    // Create SegmentQuickEdit controls
    o = new Fl_Group(0, 270, 400, 150, "SegmentQuickEdit");
    o->user_data((void*)(this));
    o->labeltype(FL_NO_LABEL);
    o->box(FL_ENGRAVED_BOX);
    {
      {
        Fl_Input *o = editSegmentLabel = new Fl_Input(10, 290, 370, 30, "Label");
        o->align(133);
        o->callback((Fl_Callback*)cb_editSegmentLabel);
      }
      {
        Fl_Slider *o = editSegmentTransperancy = new Fl_Value_Slider(10, 345, 370, 20, "Transperancy");
        o->minimum(0);
        o->maximum(1);
        o->box(FL_EMBOSSED_BOX);
        o->type(5);
        o->callback((Fl_Callback*)cb_editSegmentTransperancy);
        o->align(FL_ALIGN_TOP_LEFT);
      }
      {
        Fl_Button *o = editSegmentColor = new Fl_Button(110, 380, 80, 30, "");
        o->callback((Fl_Callback *)cb_editSegmentPickColor);
        o->type(FL_NORMAL_BUTTON);
        o->box(FL_BORDER_BOX);
        o->label("Color chooser ");
        o->align(FL_ALIGN_LEFT);
        o->color(FL_WHITE);
      }
      {
        Fl_Check_Button *o = editSegmentVisibility = new Fl_Check_Button(275, 380, 100, 30, "Display label");
        o->callback((Fl_Callback *)cb_editSegmentVisibility);
        o->align(FL_ALIGN_WRAP);
      }
    }
    o->end(); // End of SegmentQuickEdit controls

    // Creates histgram tools controls
    _histogramWindow = NULL;

    // Creates segmentation tools controls
    o = new Fl_Group(0, 420, 400, 240, "Segmentation_tools");
    o->user_data((void*)(this));
    o->labeltype(FL_NO_LABEL);
    o->box(FL_ENGRAVED_BOX);
    {
      regionGrowingMin = new Fl_Value_Slider2(10, 570, 380, 20, "Minimum threshold for region growing");
      regionGrowingMin->type(5);
      regionGrowingMin->box(FL_EMBOSSED_BOX);
      regionGrowingMin->step(1);
      regionGrowingMin->callback((Fl_Callback*)cb_SetRegionGrowingThresholdMinimum);
      regionGrowingMin->align(FL_ALIGN_BOTTOM_LEFT);
      regionGrowingMax = new Fl_Value_Slider2(10, 610, 380, 20, "Maximum threshold for region growing");
      regionGrowingMax->type(5);
      regionGrowingMax->box(FL_EMBOSSED_BOX);
      regionGrowingMax->step(1);
      regionGrowingMax->callback((Fl_Callback*)cb_SetRegionGrowingThresholdMaximum);
      regionGrowingMax->align(FL_ALIGN_BOTTOM_LEFT);
      paintBrushWidth = new Fl_Choice(49, 510, 52, 30, "Brush width");
      paintBrushWidth->align(FL_ALIGN_BOTTOM);
      paintBrushWidth->menu(menu_PaintBrushWidth);
      regionGrowingMode = new Fl_Choice(169, 510, 52, 30, "Region growing");
      regionGrowingMode->align(FL_ALIGN_BOTTOM);
      regionGrowingMode->menu(menu_RegionGrowing);
      Fl_Button* o1 = new Fl_Button(289, 440, 52, 52);
      o1->callback((Fl_Callback*)cb_showHistogram);
      Fl_Pixmap *p1 = new Fl_Pixmap(kchart_xpm);
      p1->label(o1);
      Fl_Button* o2 = drawContour = new Fl_Button(49, 440, 52, 52);
      o2->callback((Fl_Callback*)cb_DrawContour);
      Fl_Pixmap *p2 = new Fl_Pixmap(color_line_open_xpm);
      p2->label(o2);
      Fl_Button* o3 = drawPaintBrush = new Fl_Button(169, 440, 52, 52);
      o3->callback((Fl_Callback*)cb_DrawPaintBrush);
      Fl_Pixmap *p3 = new Fl_Pixmap(color_line_closed_xpm);
      p3->label(o3);
    }
    o->end();
  }
}

void Fl_RViewUI::DisplayLabels()
{
  displayLabels->value(1);
}

#endif
