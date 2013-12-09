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

#include <Fl_RViewUI.h>

Fl_RViewUI  *rviewUI;
Fl_RView    *viewer;
irtkRView    *rview;

void usage()
{
  cerr << "Usage: rview [target] <source <dofin>> <options>\n";
  cerr << "Where <options> can be one or more of the following:\n";
  cerr << "\t<-config           file.cnf>     Rview configuration file\n";
  cerr << "\t<-target_landmarks file.vtk>     Target landmarks (vtkPolyData)\n";
  cerr << "\t<-source_landmarks file.vtk>     Source landmarks (vtkPolyData)\n";
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
  cerr << "\t   where color is  <red | blue | green | rainbow>\n";
  cerr << "\t<-line value>                    Line thickness\n";  
  cerr << "\t<-diff>                          Subtraction view\n";
  cerr << "\t<-tcontour r g b>                Switch on target contours (see -tmin)\n";
  cerr << "\t<-scontour r g b>                Switch on source contours (see -smin)\n";
  cerr << "\t<-seg              file.nii.gz>  Labelled segmentation image\n";
  cerr << "\t<-lut              file.seg>     Colour lookup table for labelled segmentation\n\n";
  cerr << "\t<-labels>                        Display segmentation labels instead of contours\n";
  cerr << "\t<-background r g b>              RGB color for background (unsigned char)\n";

  
  cerr << "\tMouse events:\n";
  cerr << "\tLeft mouse click :               Reslice\n\n";
  rview->cb_keyboard_info();
  // exit(1);
}

int main(int argc, char **argv)
{
  int ok;
  int filename_argc;
  char **filename_argv;

  rviewUI = new Fl_RViewUI;

  if (argc == 1) {
    usage();
  }
  if ((argc > 1) && (argv[1][0] != '-' )) {
    rview->ReadTarget(argv[1]);
    rviewUI->AddTarget(argv[1]);
    argc--;
    argv++;
  }
  if ((argc > 1) && (argv[1][0] != '-' )) {
    rview->ReadSource(argv[1]);
    rviewUI->AddSource(argv[1]);
    argc--;
    argv++;
  }
  if ((argc > 1) && (argv[1][0] != '-' )) {
    rview->ReadTransformation(argv[1]);
    rviewUI->AddTransformation(argv[1]);
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
    if ( (ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      rview->ReadTransformation(argv[1]);
      rviewUI->AddTransformation(argv[1]);
      argc--;
      argv++;
      while ((argc > 1) && (argv[1][0] != '-' )) {
        rviewUI->AddTransformation(argv[1]);
        argc--;
        argv++;
      }
      ok = true;
    }
    if ( (ok == false) && (strcmp(argv[1], "-config") == 0)) {
      argc--;
      argv++;
      rview->Read(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ( (ok == false) && (strcmp(argv[1], "-target_landmarks") == 0)) {
      argc--;
      argv++;
      rview->ReadTargetLandmarks(argv[1]);
      rview->DisplayLandmarksOn();
      argc--;
      argv++;
      ok = true;
    }
    if ( (ok == false) && (strcmp(argv[1], "-source_landmarks") == 0)) {
      argc--;
      argv++;
      rview->ReadSourceLandmarks(argv[1]);
      rview->DisplayLandmarksOn();
      argc--;
      argv++;
      ok = true;
    }
#ifdef HAS_VTK
    if ( (ok == false) && (strcmp(argv[1], "-object") == 0)) {
      argc--;
      argv++;
      rview->ReadObject(argv[1]);
      rview->DisplayObjectOn();
      argc--;
      argv++;
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
    if ((ok == false) && (strcmp(argv[1], "-cursor") == 0)) {
      argc--;
      argv++;
      rview->DisplayCursorOff();
      rview->DisplayAxisLabelsOff();
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
    if ((ok == false) && (strcmp(argv[1], "-cspline") == 0)) {
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
      argv++;
      argc--;
      rview->DisplayTargetContoursOn();
      rview->SetTargetContourColor(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
      argv++;
      argc--;
      argv++;
      argc--;
      argv++;
      argc--;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-scontour") == 0)){
      argv++;
      argc--;
      rview->DisplaySourceContoursOn();
      rview->SetSourceContourColor(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
      argv++;
      argc--;
      argv++;
      argc--;
      argv++;
      argc--;
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-seg") == 0)){
      argv++;
      argc--;
      rview->ReadSegmentation(argv[1]);
      rviewUI->UpdateSegmentationBrowser();
      argv++;
      argc--;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-lut") == 0)){
      argv++;
      argc--;
      rview->GetSegmentTable()->Read(argv[1]);
      rview->SegmentationContoursOn();
      rviewUI->UpdateSegmentationBrowser();
      argv++;
      argc--;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-line") == 0)){
      argv++;
      argc--;
      double val = atof(argv[1]);
      if ((val > 0.9) && (val < 10.1)){
        rview->SetLineThickness(val);
      }
      argv++;
      argc--;
      ok = true;
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
    if ((ok == false) && (strcmp(argv[1], "-labels") == 0)){
      argv++;
      argc--;
      rview->SegmentationContoursOff();
      rview->SegmentationLabelsOn();
      rview->SegmentationUpdateOn();
      ok = true;
    }    
    // Ignore the following display specific options
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
      argv++;
      argc--;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-background") == 0)) {
      argv++;
      argc--;
      rview->SetBackgroundColor(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));
      argv++;
      argc--;
      argv++;
      argc--;
      argv++;
      argc--;      
      ok = true;
    }
    
    if (ok == false) {
      cerr << "Unknown argument: " << argv[1] << endl;
      exit(1);
    }
  }

  Fl::visual(FL_DOUBLE | FL_RGB);
  rviewUI->show();
  return Fl::run();
}
