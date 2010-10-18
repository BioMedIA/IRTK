/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkTransformation.h>

#include <irtkRegistration.h>

#ifdef __APPLE__
#include <OpenGl/gl.h>
#include <OpenGl/glu.h>
#include <glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <irtkRView.h>

irtkRView *rview;

char *offscreen_file;
int target_min, target_max, target_delta;
int source_min, source_max, source_delta;

void usage ()
{
  cerr << "Usage: display [target] <source <dofin>> <options>\n";
  cerr << "Where <options> can be one or more of the following:\n";
  cerr << "\t<-config           file.cnf>     Rview configuration file\n";
  cerr << "\t<-target_landmarks file.vtk>     Target Landmarks (vtkPolyData)\n";
  cerr << "\t<-source_landmarks file.vtk>     Source Landmarks (vtkPolyData)\n";
  cerr << "\t<-target_isolines>               Target isolines\n";
  cerr << "\t<-source_isolines>               Source isolines\n";
#ifdef HAS_VTK
  cerr << "\t<-object           file.vtk>     Object           (vtkPointSet)\n";
  cerr << "\t<-object_warp>                   Warp object with vectors\n";
  cerr << "\t<-object_grid>                   Object grid on\n";
#endif
  cerr << "\t<-eigen values.irtk vectors.irtk>  Eigen modes\n";
  cerr << "\t<-xy      | -xz      | -yz>      Single     view\n";
  cerr << "\t<-xy_xz_v | -xy_yz_v | -xz_yz_v> Vertical   view\n";
  cerr << "\t<-xy_xz_h | -xy_yz_h | -xz_yz_h> Horizontal view\n";
  cerr << "\t<-cursor>                        Cursor off\n";
  cerr << "\t<-grid>                          Deformation grid   on\n";
  cerr << "\t<-points>                        Deformation points on\n";
  cerr << "\t<-arrow>                         Deformation arrows on\n";
  cerr << "\t<-level value>                   Deformation level\n";
  cerr << "\t<-res   value>                   Resolution factor\n";
  cerr << "\t<-nn>                            Nearest neighbour interpolation (default)\n";
  cerr << "\t<-linear>                        Linear interpolation\n";
  cerr << "\t<-c1spline>                      C1-spline interpolation\n";
  cerr << "\t<-bspline>                       B-spline interpolation\n";
  cerr << "\t<-sinc>                          Sinc interpolation\n";
  cerr << "\t<-origin x y z>                  Reslice position\n";
  cerr << "\t<-tmin value>                    Min. target intensity\n";
  cerr << "\t<-tmax value>                    Max. target intensity\n";
  cerr << "\t<-smin value>                    Min. source intensity\n";
  cerr << "\t<-smax value>                    Max. source intensity\n";
  cerr << "\t<-sub_min value>                 Min. subtraction intensity\n";
  cerr << "\t<-sub_max value>                 Max. subtraction intensity\n";
  cerr << "\t<-view_target>                   View target (default)\n";
  cerr << "\t<-view_source>                   View source\n";
  cerr << "\t<-mix>                           Mixed viewport (checkerboard)\n";
  cerr << "\t<-tcolor color>                  Target image color\n";
  cerr << "\t<-scolor color>                  Source image color\n";
  cerr << "\t<-line value>                    Line thickness\n";
  cerr << "\t   where color is  <red | blue | green | rainbow>\n";
  cerr << "\t<-diff>                          Subtraction view\n";
  cerr << "\t<-tcontour>                      Switch on target contours (see -tmin)\n";
  cerr << "\t<-scontour>                      Switch on source contours (see -smin)\n";
  cerr << "\t<-seg              file.nii.gz>  Labelled segmentation image\n";
  cerr << "\t<-lut              file.seg>     Colour lookup table for labelled segmentation\n\n";
  cerr << "\tDisplay specific options:\n";
  cerr << "\t<-x value>                       Width\n";
  cerr << "\t<-y value>                       Height\n";
  cerr << "\t<-offscreen>                     Offscreen rendering\n\n";
  cerr << "\tMouse events:\n";
  cerr << "\tLeft mouse click :               Reslice\n\n";
  rview->cb_special_info();
  rview->cb_keyboard_info();
  exit(1);
}

void display(void)
{
  rview->Draw();
  glutSwapBuffers();
}

void reshape (int w, int h)
{
  rview->Resize(w, h);
}

void mouse(int button, int state, int x, int y)
{
  double a, b, c;

  if (state == GLUT_UP) return;

  if (button == GLUT_LEFT_BUTTON) {
    rview->SetOrigin(x, y);
    rview->Update();
    rview->Draw();
    glutSwapBuffers();

    // Say what we have done
    rview->GetOrigin(a, b, c);
    cout << "Resliced at " << a << " " << b << " " << c << endl;
  }
}

void special(int key, int x, int y)
{
  rview->cb_special(key, x, y, target_delta, source_delta);
  glutSwapBuffers();

  return;
}

void keyboard(unsigned char key, int, int)
{
  rview->cb_keyboard(key);
  glutSwapBuffers();

  return;
}

int main(int argc, char** argv)
{
  int i, x, y;
  bool ok, offscreen;
  int filename_argc;
  char **filename_argv;

  x = 512;
  y = 512;
  offscreen = false;
  for (i = 0; i < argc-1; i++) {
    if (strcmp(argv[i], "-x") == 0) {
      x = atoi(argv[i+1]);
      i++;
    }
    if (strcmp(argv[i], "-y") == 0) {
      y = atoi(argv[i+1]);
      i++;
    }
  }

  // Initialize viewer
  rview = new irtkRView(x, y);

  if (argc == 1) {
    usage();
  }
  if ((argc > 1) && (argv[1][0] != '-' )) {
    rview->ReadTarget(argv[1]);
    argc--;
    argv++;
  }
  if ((argc > 1) && (argv[1][0] != '-' )) {
    rview->ReadSource(argv[1]);
    argc--;
    argv++;
  }
  if ((argc > 1) && (argv[1][0] != '-' )) {
    rview->ReadTransformation(argv[1]);
    argc--;
    argv++;
  }

  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-target") == 0)) {
      argc--;
      argv++;

      // Get filenames for sequence
      filename_argc = 0;
      filename_argv = argv;
      while ((argc > 1) && (argv[1][0] != '-' )) {
        argv++;
        argc--;
        filename_argc++;
      }
      filename_argv++;

      // Read sequence
      rview->ReadTarget(filename_argc, filename_argv);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-source") == 0)) {
      argc--;
      argv++;

      // Get filenames for sequence
      filename_argc = 0;
      filename_argv = argv;
      while ((argc > 1) && (argv[1][0] != '-' )) {
        argv++;
        argc--;
        filename_argc++;
      }
      filename_argv++;

      // Read sequence
      rview->ReadSource(filename_argc, filename_argv);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      rview->ReadTransformation(argv[1]);
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-config") == 0)) {
      argc--;
      argv++;
      rview->Read(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-target_landmarks") == 0)) {
      argc--;
      argv++;
      rview->ReadTargetLandmarks(argv[1]);
      rview->DisplayLandmarksOn();
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-source_landmarks") == 0)) {
      argc--;
      argv++;
      rview->ReadSourceLandmarks(argv[1]);
      rview->DisplayLandmarksOn();
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-target_isolines") == 0)) {
      argc--;
      argv++;
      rview->DisplayTargetContoursOn();
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-source_isolines") == 0)) {
      argc--;
      argv++;
      rview->DisplaySourceContoursOn();
      ok = true;
    }
#ifdef HAS_VTK
    if ((ok == false) && (strcmp(argv[1], "-object") == 0)) {
      argc--;
      argv++;
      rview->ReadObject(argv[1]);
      rview->DisplayObjectOn();
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-object_warp") == 0)) {
      argc--;
      argv++;
      rview->DisplayObjectWarpOn();
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-object_grid") == 0)) {
      argc--;
      argv++;
      rview->DisplayObjectGridOn();
      ok = true;
    }
#endif
    if ((ok == false) && (strcmp(argv[1], "-xy") == 0)) {
      argc--;
      argv++;
      rview->Configure(View_XY);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xz") == 0)) {
      argc--;
      argv++;
      rview->Configure(View_XZ);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-yz") == 0)) {
      argc--;
      argv++;
      rview->Configure(View_YZ);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xy_xz_v") == 0)) {
      argc--;
      argv++;
      rview->Configure(View_XY_XZ_v);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xy_yz_v") == 0)) {
      argc--;
      argv++;
      rview->Configure(View_XY_YZ_v);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xz_yz_v") == 0)) {
      argc--;
      argv++;
      rview->Configure(View_XZ_YZ_v);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xy_xz_h") == 0)) {
      argc--;
      argv++;
      rview->Configure(View_XY_XZ_h);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xy_yz_h") == 0)) {
      argc--;
      argv++;
      rview->Configure(View_XY_YZ_h);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xz_yz_h") == 0)) {
      argc--;
      argv++;
      rview->Configure(View_XZ_YZ_h);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-labels") == 0)) {
      argc--;
      argv++;
      rview->DisplayAxisLabelsOff();
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-cursor") == 0)) {
      argc--;
      argv++;
      rview->DisplayCursorOff();
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-grid") == 0)) {
      argc--;
      argv++;
      rview->DisplayDeformationGridOn();
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-points") == 0)) {
      argc--;
      argv++;
      rview->DisplayDeformationPointsOn();
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-arrow") == 0)) {
      argc--;
      argv++;
      rview->DisplayDeformationArrowsOn();
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-tmax") == 0)) {
      argc--;
      argv++;
      rview->SetDisplayMaxTarget(atof(argv[1]));
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-tmin") == 0)) {
      argc--;
      argv++;
      rview->SetDisplayMinTarget(atof(argv[1]));
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-smax") == 0)) {
      argc--;
      argv++;
      rview->SetDisplayMaxSource(atof(argv[1]));
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-smin") == 0)) {
      argc--;
      argv++;
      rview->SetDisplayMinSource(atof(argv[1]));
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sub_max") == 0)) {
      argc--;
      argv++;
      rview->SetDisplayMaxSubtraction(atof(argv[1]));
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sub_min") == 0)) {
      argc--;
      argv++;
      rview->SetDisplayMinSubtraction(atof(argv[1]));
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-res") == 0)) {
      argc--;
      argv++;
      rview->SetResolution(atof(argv[1]));
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-origin") == 0)) {
      argc--;
      argv++;
      rview->SetOrigin(atof(argv[1]), atof(argv[2]), atof(argv[3]));
      argc--;
      argv++;
      argc--;
      argv++;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-nn") == 0)) {
      rview->SetTargetInterpolationMode(Interpolation_NN);
      rview->SetSourceInterpolationMode(Interpolation_NN);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-linear") == 0)) {
      rview->SetTargetInterpolationMode(Interpolation_Linear);
      rview->SetSourceInterpolationMode(Interpolation_Linear);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-c1spline") == 0)) {
      rview->SetTargetInterpolationMode(Interpolation_CSpline);
      rview->SetSourceInterpolationMode(Interpolation_CSpline);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-bspline") == 0)) {
      rview->SetTargetInterpolationMode(Interpolation_BSpline);
      rview->SetSourceInterpolationMode(Interpolation_BSpline);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-sinc") == 0)) {
      rview->SetTargetInterpolationMode(Interpolation_Sinc);
      rview->SetSourceInterpolationMode(Interpolation_Sinc);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-view_target") == 0)) {
      rview->SetViewMode(View_A);
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-view_source") == 0)) {
      rview->SetViewMode(View_B);
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mix") == 0)) {
      rview->SetViewMode(View_Checkerboard);
      argv++;
      argc--;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-diff") == 0)){
      rview->SetViewMode(View_Subtraction);
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-tcontour") == 0)){
      rview->DisplayTargetContoursOn();
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-scontour") == 0)){
      rview->DisplaySourceContoursOn();
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-seg") == 0)){
      argv++;
      argc--;
      rview->ReadSegmentation(argv[1]);
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-lut") == 0)){
      argv++;
      argc--;
      rview->GetSegmentTable()->Read(argv[1]);
      rview->SegmentationContoursOn();
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-labels") == 0)){
      argv++;
      argc--;
      rview->SegmentationContoursOff();
      rview->SegmentationLabelsOn();
      rview->SegmentationUpdateOn();
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-line") == 0)) {
      argc--;
      argv++;
	  rview->SetLineThickness(atoi(argv[1]));
	  argc--;
	  argv++;
	}
	if ((ok == false) && (strcmp(argv[1], "-tcolor") == 0)) {
      argc--;
      argv++;
      if (strcmp(argv[1], "red") == 0) {
        rview->GetTargetLookupTable()->SetColorModeToRed();
        argv++;
        argc--;
        ok = true;
      } else {
        if (strcmp(argv[1], "green") == 0) {
          rview->GetTargetLookupTable()->SetColorModeToGreen();
          argv++;
          argc--;
          ok = true;
        } else {
          if (strcmp(argv[1], "blue") == 0) {
            rview->GetTargetLookupTable()->SetColorModeToBlue();
            argv++;
            argc--;
            ok = true;
          } else {
            if (strcmp(argv[1], "rainbow") == 0) {
              rview->GetTargetLookupTable()->SetColorModeToRainbow();
              argv++;
              argc--;
              ok = true;
            }
          }
        }
      }
    }
    if ((ok == false) && (strcmp(argv[1], "-scolor") == 0)) {
      argc--;
      argv++;
      if (strcmp(argv[1], "red") == 0) {
        rview->GetSourceLookupTable()->SetColorModeToRed();
        argv++;
        argc--;
        ok = true;
      } else {
        if (strcmp(argv[1], "green") == 0) {
          rview->GetSourceLookupTable()->SetColorModeToGreen();
          argv++;
          argc--;
          ok = true;
        } else {
          if (strcmp(argv[1], "blue") == 0) {
            rview->GetSourceLookupTable()->SetColorModeToBlue();
            argv++;
            argc--;
            ok = true;
          } else {
            if (strcmp(argv[1], "rainbow") == 0) {
              rview->GetSourceLookupTable()->SetColorModeToRainbow();
              argv++;
              argc--;
              ok = true;
            }
          }
        }
      }
    }
    // Display specific options
    if ((ok == false) && (strcmp(argv[1], "-x") == 0)) {
      argc--;
      argv++;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-y") == 0)) {
      argc--;
      argv++;
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-offscreen") == 0)) {
      offscreen = true;
      argv++;
      argc--;
      offscreen_file = argv[1];
      argv++;
      argc--;
      ok = true;
    }
    if (ok == false) {
      cerr << "Unknown argument: " << argv[1] << endl;
      exit(1);
    }
  }

  // Initilaize min/max/delta values for special function keys
  target_min   = rview->GetTargetMin();
  target_max   = rview->GetTargetMax();
  target_delta = round((target_max - target_min) / 50.0);
  source_min   = rview->GetSourceMin();
  source_max   = rview->GetSourceMax();
  source_delta = round((source_max - source_min) / 50.0);

  // Initialize graphics window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(rview->GetWidth(), rview->GetHeight());
  glutInitWindowPosition(0, 0);
  glutCreateWindow("RView");

  if (offscreen == true) {

    // Start rendering to file
    rview->DrawOffscreen(offscreen_file);

  } else {

    // Initialize callback functions
    glutMouseFunc(mouse);
    glutSpecialFunc(special);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);

    // Start rendering
    glutMainLoop();
  }

  return 0;
}
