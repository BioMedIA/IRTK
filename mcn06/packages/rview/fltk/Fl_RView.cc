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

#include <Fl_RView.h>

#include <Fl_RViewUI.h>

extern Fl_RViewUI *rviewUI;

Fl_RView::Fl_RView(int x, int y, int w, int h, const char *name) : Fl_Gl_Window(x, y, w, h, name)
{
  v = new irtkRView(w, h);
}

void Fl_RView::draw()
{
  if (!valid()) {
    v->Resize(w(), h());
  }
  v->Draw();
}

int Fl_RView::handle(int event)
{
  char buffer1[256], buffer2[256], buffer3[256], buffer4[256], buffer5[256];

  switch (event) {
  case FL_KEYBOARD:
  case FL_SHORTCUT:
    if (Fl::event_key() == FL_Delete) {
      v->ClearContour();
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if ((Fl::event_key() == FL_BackSpace) && (Fl::event_shift() != 0)) {
      v->ClearContour();
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if (Fl::event_key() == FL_BackSpace) {
      v->UndoContour();
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if (Fl::event_key() == FL_Enter) {
#ifdef HAS_SEGMENTATION_PANEL
      if (rviewUI->_selected == 0) {
        fl_alert("Please, select a label!\n");
        return 1;
      } else {
        v->FillContour(rviewUI->_id, 0);
        v->ClearContour();
        v->SegmentationLabelsOn();
        v->SegmentationUpdateOn();
        v->Update();
#ifdef HAS_VISUALISATION_PANEL
        rviewUI->UpdateImageControlWindow();
#endif
        rviewUI->update();
        this->redraw();
        return 1;
      }
#endif
    }
    if (Fl::event_key() == FL_Tab) {
#ifdef HAS_SEGMENTATION_PANEL
      v->RegionGrowContour(Fl::event_x(), Fl::event_y());
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
#endif
    }
    v->cb_keyboard(Fl::event_text()[0]);
    rviewUI->update();
    this->redraw();
    return 1;
    break;
  case FL_PUSH:
    if ((Fl::event_button() == 1) && (Fl::event_shift() != 0)) {
      v->AddContour(Fl::event_x(), Fl::event_y(), FirstPoint);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
#ifdef HAS_SEGMENTATION_PANEL
    if ((Fl::event_button() == 3) && (Fl::event_shift() != 0)) {
      v->FillArea(Fl::event_x(), Fl::event_y());
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
#endif
    if ((Fl::event_button() == 1) && (Fl::event_ctrl() != 0)) {
      v->UpdateROI1(Fl::event_x(), Fl::event_y());
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if ((Fl::event_button() == 3) && (Fl::event_ctrl() != 0)) {
      v->UpdateROI2(Fl::event_x(), Fl::event_y());
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if (Fl::event_button() == 1) {
      v->SetOrigin(Fl::event_x(), Fl::event_y());
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    break;
  case FL_DRAG:
    if ((Fl::event_state() & (FL_BUTTON1 | FL_SHIFT)) == (FL_BUTTON1 | FL_SHIFT)) {
      v->AddContour(Fl::event_x(), Fl::event_y(), NewPoint);
      v->Update();
      rviewUI->info_voxel->value(buffer1);
      rviewUI->info_world->value(buffer2);
      rviewUI->info_target->value(buffer3);
      rviewUI->info_source->value(buffer4);
      rviewUI->info_segmentation->value(buffer5);
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if ((Fl::event_state() & (FL_BUTTON1 | FL_CTRL)) == (FL_BUTTON1 | FL_CTRL)) {
      v->UpdateROI1(Fl::event_x(), Fl::event_y());
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    if ((Fl::event_state() & (FL_BUTTON3 | FL_CTRL)) == (FL_BUTTON3 | FL_CTRL)) {
      v->UpdateROI2(Fl::event_x(), Fl::event_y());
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    break;
  case FL_RELEASE:
    if ((Fl::event_button() == 1) && (Fl::event_shift() != 0)) {
      v->AddContour(Fl::event_x(), Fl::event_y(), LastPoint);
      v->Update();
      rviewUI->update();
      this->redraw();
      return 1;
    }
    return 0;
  case FL_MOVE:
    v->MousePosition(Fl::event_x(), Fl::event_y());
    rviewUI->update();
    this->redraw();
    return 1;
    break;
  case FL_MOUSEWHEEL:
    v->MouseWheel(Fl::event_x(), Fl::event_y(), Fl::event_dy());
    v->Update();
    rviewUI->update();
    this->redraw();
    return 1;
    break;
  case FL_FOCUS:
    return 1;
    break;
  case FL_UNFOCUS:
    return 1;
    break;
  case FL_ENTER:
    Fl::focus(this);
    return 1;
    break;
  case FL_LEAVE:
    return 1;
    break;
  default:
    return Fl_Gl_Window::handle(event);
  }
  return 0;
}

