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

#ifdef HAS_REGISTRATION_PANEL

extern Fl_RViewUI  *rviewUI;
extern Fl_RView    *viewer;
extern irtkRView    *rview;

void registration_cb1()
{
  // Wait for any interaction
  Fl::wait(0);

  // Update display
  rview->SourceUpdateOn();
  rview->Update();

  // Update
  viewer->redraw();
  rviewUI->update();
}

void registration_cb2()
{
  // Wait for any interaction
  Fl::wait(0);

  // Update
  viewer->redraw();
  rviewUI->update();
}

Fl_Menu_Item Fl_RViewUI::menu_similarityMeasure[] = {
  {"JE", 0, (Fl_Callback*)cb_similarityMeasure, "JE"},
  {"MI", 0, (Fl_Callback*)cb_similarityMeasure, "MI"},
  {"NMI", 0, (Fl_Callback*)cb_similarityMeasure, "NMI"},
  {"CR_XY", 0, (Fl_Callback*)cb_similarityMeasure, "CR_XY"},
  {"CR_YX", 0, (Fl_Callback*)cb_similarityMeasure, "CR_YX"},
  {"CC", 0, (Fl_Callback*)cb_similarityMeasure, "CC"},
  {"SSD", 0, (Fl_Callback*)cb_similarityMeasure, "SSD"},
  {0}
};

void Fl_RViewUI::cb_similarityMeasure(Fl_Widget *, char *v)
{
  if (strcmp(v, "JE") == 0) {
    rview->GetRegistration()->PutSimilarityMeasure(JE);
  }
  if (strcmp(v, "MI") == 0) {
    rview->GetRegistration()->PutSimilarityMeasure(MI);
  }
  if (strcmp(v, "NMI") == 0) {
    rview->GetRegistration()->PutSimilarityMeasure(NMI);
  }
  if (strcmp(v, "CR_XY") == 0) {
    rview->GetRegistration()->PutSimilarityMeasure(CR_XY);
  }
  if (strcmp(v, "CR_YX") == 0) {
    rview->GetRegistration()->PutSimilarityMeasure(CR_YX);
  }
  if (strcmp(v, "CC") == 0) {
    rview->GetRegistration()->PutSimilarityMeasure(CC);
  }
  if (strcmp(v, "SSD") == 0) {
    rview->GetRegistration()->PutSimilarityMeasure(SSD);
  }
}

void Fl_RViewUI::cb_numberOfBins(Fl_Counter *o, void *)
{
  rview->GetRegistration()->PutNumberOfBins(round(o->value()));
}

void Fl_RViewUI::cb_numberOfIterations(Fl_Counter *o, void *)
{
  rview->GetRegistration()->PutNumberOfIterations(round(o->value()));
}

void Fl_RViewUI::cb_numberOfSteps(Fl_Counter *o, void *)
{
  rview->GetRegistration()->PutNumberOfSteps(round(o->value()));
}

void Fl_RViewUI::cb_lengthOfSteps(Fl_Counter *o, void *)
{
  rview->GetRegistration()->PutLengthOfSteps(o->value());
}

void Fl_RViewUI::cb_numberOfLevels(Fl_Button *, char *v)
{
  rview->GetRegistration()->PutNumberOfLevels(atoi(v));
}

void Fl_RViewUI::cb_targetBlurring(Fl_Value_Slider *, void *)
{
  rview->GetRegistration()->PutTargetBlurring(rviewUI->targetBlurring->value());
}

void Fl_RViewUI::cb_targetResolution(Fl_Value_Slider *, void *)
{
  rview->GetRegistration()->PutTargetResolution(rviewUI->targetResolution->value());
}

void Fl_RViewUI::cb_sourceBlurring(Fl_Value_Slider *, void *)
{
  rview->GetRegistration()->PutSourceBlurring(rviewUI->sourceBlurring->value());
}

void Fl_RViewUI::cb_sourceResolution(Fl_Value_Slider *, void *)
{
  rview->GetRegistration()->PutSourceResolution(rviewUI->sourceResolution->value());
}

void Fl_RViewUI::cb_targetPadding(Fl_Value_Slider *, void *)
{
  rview->GetRegistration()->PutTargetPadding(round(rviewUI->targetPadding->value()-1));
}

void Fl_RViewUI::cb_loadParameter(Fl_Menu_Bar *, void *v)
{
  char *filename = fl_file_chooser("Load parameter", "*.*", "");

  if (filename != NULL) {
    rview->ReadParameter(filename);
    rviewUI->update();
  }
}

void Fl_RViewUI::cb_saveParameter(Fl_Menu_Bar *, void *v)
{
  char *filename = fl_file_chooser("Save parameter", "*.*", "");

  if (filename != NULL) {
    rview->WriteParameter(filename);
  }
}

void Fl_RViewUI::cb_startRegistration(Fl_Button* o, void* v)
{
  int i;
  double x1, y1, z1, x2, y2, z2;
  char buffer[256];

  if ((rview->GetTarget()->IsEmpty() != True) &&
      (rview->GetSource()->IsEmpty() != True)) {

    irtkRegistration *registration = rview->GetRegistration();

    // Calculate ROI
    rview->GetROI(x1, y1, z1, x2, y2, z2);
    rview->GetTarget()->WorldToImage(x1, y1, z1);
    rview->GetTarget()->WorldToImage(x2, y2, z2);

    // Copy images
    irtkGreyImage target = rview->GetTarget()->GetRegion(round(x1), round(y1), round(z1), round(x2)+1, round(y2)+1, round(z2)+1);
    irtkGreyImage source = *rview->GetSource();

    // Run registration
    registration->SetInput(&target, &source);
    registration->SetOutput(rview->GetTransformation());
    registration->SetCallback1(registration_cb1);
    registration->SetCallback2(registration_cb2);
    registration->Run();

    // Update transformation browser in case things have changed
    for (i = 0; i <= rview->NumberOfDeformationLevels(); i++) {
      rview->GetTransformationText(i, buffer);
      rviewUI->transformationBrowser->text(i+1, buffer);
    }

    // Alert user
    fl_alert("Registration has finished");

    // Delete registration
    delete registration;
  }
}

void Fl_RViewUI::cb_stopRegistration(Fl_Button* o, void* v)
{
  // Alert user
  fl_alert("You can't stop the registration");
}

void Fl_RViewUI::cb_guessRegistrationParameter(Fl_Button* o, void* v)
{
  rview->GuessRegistrationParameter();
  rviewUI->UpdateRegistrationControlWindow();

  // Update
  viewer->redraw();
  rviewUI->update();
}

void Fl_RViewUI::UpdateRegistrationControlWindow()
{
  int i;

  targetBlurring->value(rview->GetRegistration()->GetTargetBlurring());
  targetResolution->value(rview->GetRegistration()->GetTargetResolution());
  sourceBlurring->value(rview->GetRegistration()->GetSourceBlurring());
  sourceResolution->value(rview->GetRegistration()->GetSourceResolution());
  numberOfBins->value(rview->GetRegistration()->GetNumberOfBins());
  numberOfIterations->value(rview->GetRegistration()->GetNumberOfIterations());
  numberOfSteps->value(rview->GetRegistration()->GetNumberOfSteps());
  lengthOfSteps->value(rview->GetRegistration()->GetLengthOfSteps());
  targetPadding->minimum(rview->GetTargetLookupTable()->minData);
  targetPadding->maximum(rview->GetTargetLookupTable()->maxData);
  if (rview->GetRegistration()->GetTargetPadding() <
      targetPadding->minimum()) {
    targetPadding->value(targetPadding->minimum());
  } else if (rview->GetRegistration()->GetTargetPadding() >
             targetPadding->maximum()) {
    targetPadding->value(targetPadding->maximum());
  } else {
    // By default, this is MIN_GREY
    targetPadding->value(rview->GetRegistration()->GetTargetPadding());
  }

  for (i = 0; i < 5; i++) level[i]->value(0);
  switch (rview->GetRegistration()->GetNumberOfLevels()) {
  case 0:
    level[0]->value(1);
    break;
  case 1:
    level[0]->value(1);
    break;
  case 2:
    level[1]->value(1);
    break;
  case 3:
    level[2]->value(1);
    break;
  case 4:
    level[3]->value(1);
    break;
  case 5:
    level[4]->value(1);
    break;
  default:
    break;
  }

  switch (rview->GetRegistration()->GetSimilarityMeasure()) {
  case JE:
    similarityMeasure->value(0);
    break;
  case MI:
    similarityMeasure->value(1);
    break;
  case NMI:
    similarityMeasure->value(2);
    break;
  case CR_XY:
    similarityMeasure->value(3);
    break;
  case CR_YX:
    similarityMeasure->value(4);
    break;
  case CC:
    similarityMeasure->value(5);
    break;
  case SSD:
    similarityMeasure->value(6);
    break;
  default:
    cerr << "Unknown similarity measure" << endl;
    exit(1);
  }
}

void Fl_RViewUI::InitializeRegistrationControlWindow()
{
  {
    // Create registration optimization parameters
    Fl_Window* o = new Fl_Window(5, 30, 310, 170);
    o->box(FL_ENGRAVED_BOX);
    {
      Fl_Counter* o = numberOfBins =
                        new Fl_Counter(120, 10, 180, 30, "No. of bins");
      o->callback((Fl_Callback*)cb_numberOfBins);
      o->align(FL_ALIGN_LEFT);
      o->minimum(0);
      o->maximum(1000);
      o->step(1);
      o->lstep(10);
    }
    {
      Fl_Counter* o = numberOfIterations =
                        new Fl_Counter(120, 50, 180, 30, "No. of iterations");
      o->callback((Fl_Callback*)cb_numberOfIterations);
      o->align(FL_ALIGN_LEFT);
      o->minimum(0);
      o->maximum(1000);
      o->step(1);
      o->lstep(10);
    }
    {
      Fl_Counter* o = numberOfSteps =
                        new Fl_Counter(120, 90, 180, 30, "No. of steps");
      o->callback((Fl_Callback*)cb_numberOfSteps);
      o->align(FL_ALIGN_LEFT);
      o->minimum(0);
      o->maximum(100);
      o->step(1);
      o->lstep(5);
    }
    {
      Fl_Counter* o = lengthOfSteps =
                        new Fl_Counter(120, 130, 180, 30, "Length of steps");
      o->callback((Fl_Callback*)cb_lengthOfSteps);
      o->minimum(0);
      o->align(FL_ALIGN_LEFT);
      o->maximum(20);
      o->step(0.1);
      o->lstep(1);
    }
    o->end(); // End of registration optimization parameters
  }
  {
    // Create registration multiresolution and similarity parameters
    Fl_Window* o = new Fl_Window(5, 200, 310, 80);
    o->box(FL_ENGRAVED_BOX);
    {
      Fl_Box* o;
      o = new Fl_Box(10, 10, 100, 30, "Levels");
    }
    {
      Fl_Check_Button* o = level[0] =
                             new Fl_Check_Button(100, 10, 40, 30, "1");
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->type(FL_RADIO_BUTTON);
      o->callback((Fl_Callback*)cb_numberOfLevels, "1");
    }
    {
      Fl_Check_Button* o = level[1] =
                             new Fl_Check_Button(140, 10, 40, 30, "2");
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->type(FL_RADIO_BUTTON);
      o->callback((Fl_Callback*)cb_numberOfLevels, "2");
    }
    {
      Fl_Check_Button* o = level[2] =
                             new Fl_Check_Button(180, 10, 40, 30, "3");
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->type(FL_RADIO_BUTTON);
      o->callback((Fl_Callback*)cb_numberOfLevels, "3");
    }
    {
      Fl_Check_Button* o = level[3] =
                             new Fl_Check_Button(220, 10, 40, 30, "4");
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->type(FL_RADIO_BUTTON);
      o->callback((Fl_Callback*)cb_numberOfLevels, "4");
    }
    {
      Fl_Check_Button* o = level[4] =
                             new Fl_Check_Button(260, 10, 40, 30, "5");
      o->down_box(FL_DIAMOND_DOWN_BOX);
      o->type(FL_RADIO_BUTTON);
      o->callback((Fl_Callback*)cb_numberOfLevels, "5");
    }
    {
      Fl_Choice* o = similarityMeasure =
                       new Fl_Choice(140, 45, 70, 30);
      o->align(FL_ALIGN_LEFT);
      o->label("Similarity measure: ");
      o->menu(menu_similarityMeasure);
    }
    o->end(); // End of registration multiresolution and similarity parameters
  }
  {
    // Create registration image preprocessing parameters
    Fl_Window* o = new Fl_Window(5, 280, 310, 240);
    o->box(FL_ENGRAVED_BOX);
    {
      Fl_Value_Slider* o = targetBlurring = new
      Fl_Value_Slider(10, 25, 290, 20, "Target blurring (in mm):");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->step(0.1);
      o->minimum(0);
      o->maximum(5);
      o->callback((Fl_Callback*)cb_targetBlurring);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    {
      Fl_Value_Slider* o = targetResolution = new
      Fl_Value_Slider(10, 70, 290, 20, "Target resolution (in mm):");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->step(0.1);
      o->minimum(0);
      o->maximum(5);
      o->callback((Fl_Callback*)cb_targetResolution);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    {
      Fl_Value_Slider* o = sourceBlurring = new
      Fl_Value_Slider(10, 115, 290, 20, "Source blurring (in mm):");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->step(0.1);
      o->minimum(0);
      o->maximum(5);
      o->callback((Fl_Callback*)cb_sourceBlurring);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    {
      Fl_Value_Slider* o = sourceResolution = new
      Fl_Value_Slider(10, 160, 290, 20, "Source resolution (in mm):");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->step(0.1);
      o->minimum(0);
      o->maximum(5);
      o->callback((Fl_Callback*)cb_sourceResolution);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    {
      Fl_Value_Slider* o = targetPadding = new
      Fl_Value_Slider(10, 205, 290, 20, "Target padding (ignore intensities below):");
      o->type(5);
      o->box(FL_EMBOSSED_BOX);
      o->step(1);
      o->callback((Fl_Callback*)cb_targetPadding);
      o->align(FL_ALIGN_TOP_LEFT);
    }
    o->end(); // End of registration image preprocessing parameters
  }
  {
    // Create registration controls
    Fl_Window* o = new Fl_Window(5, 520, 310, 40, "Registration controls");
    o->user_data((void*)(this));
    o->box(FL_ENGRAVED_BOX);
    {
      Fl_Button  *o = new Fl_Return_Button( 10, 5,  90,  30, "Start");
      o->callback((Fl_Callback*)cb_startRegistration);
    }
    {
      Fl_Button  *o = new Fl_Return_Button(110, 5,  90,  30, "Stop");
      o->callback((Fl_Callback*)cb_stopRegistration);
    }
    {
      Fl_Button  *o = new Fl_Return_Button(210, 5,  90,  30, "Guess");
      o->callback((Fl_Callback*)cb_guessRegistrationParameter);
    }
    o->end(); // End of registration controls
  }
}

#endif
