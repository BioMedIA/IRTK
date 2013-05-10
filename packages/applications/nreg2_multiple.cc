/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2009 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration2.h>

// Default filenames
char **source_name = NULL, **target_name = NULL;
char *dofin_name  = NULL, *dofout_name = NULL;
char *parin_name  = NULL, *parout_name = NULL;
char *mask_name = NULL;

void usage()
{
  cerr << "Usage: rreg2_multiple [Number of images] [target1...targetN] [source1...sourceN] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parin file>        Read parameter from file" << endl;
  cerr << "<-parout file>       Write parameter to file" << endl;
  cerr << "<-dofin  file>       Read transformation from file" << endl;
  cerr << "<-dofout file>       Write transformation to file" << endl;
  cerr << "<-Tp  value>         Padding value in target" << endl;
  cerr << "<-ds  value>         Initial control point spacing" << endl;
  cerr << "<-Sp  value>         Smoothness preservation value" << endl;
  cerr << "                     before running registration." << endl;
  cerr << "<-debug>             Enable debugging information" << endl;
  cerr << "<-mask file>         Use a mask to define the ROI." << endl;
  cerr << "                     Volume perservation is applied to" <<endl;
  cerr << "                     Voxels in the mask." << endl;
  cerr << "<-mask_dilation n>   Dilate mask n times before using it" << endl;
  cerr << "<-weighting levelnumber weight1...weightn>    weighting for the images" << endl;
#ifdef HAS_VTK
  cerr << "<-landmarks t s>     Landmark Regulation target and source" << endl;
#endif

  exit(1);
}

int main(int argc, char **argv)
{
  int ok, padding;
  int i, j, numberOfImages;
  int mask_dilation = false;
  double sp,spacing,**weight;

  sp = 0;
  spacing = 0;

#ifdef HAS_VTK
  vtkPolyData **landmarks;
  vtkPolyDataReader *reader;

  landmarks = NULL;
  reader = NULL;
#endif

  // Check command line
  if (argc < 3) {
    usage();
  }

  // Get image names for sequence
  numberOfImages = atoi(argv[1]);
  argv++;
  argc--;

  // Initialize the weighting
  weight = new double*[10];

  for (i = 0; i < 10; i++) {
      weight[i] = new double[numberOfImages];
      for (j = 0; j < numberOfImages; j++) {
          weight[i][j] = 1;
      }
  }

  target_name = &argv[1];
  for(i = 0; i < numberOfImages; i++){
      if(argv[1][0] == '-' ){
          usage();
          exit(1);
      }
      argv++;
      argc--;
  }
  source_name = &argv[1];
  for(i = 0; i < numberOfImages; i++){
      if(argv[1][0] == '-' ){
          usage();
          exit(1);
      }
      argv++;
      argc--;
  }

  // Read targets
  cout << "Reading targets ... "; cout.flush();
  irtkGreyImage **target = new irtkGreyImage *[numberOfImages];
  for (i = 0; i < numberOfImages; i++) {
      cout << target_name[i] << endl;
      target[i] = new irtkGreyImage(target_name[i]);
  }
  cout << "done" << endl;

  // Read sources
  cout << "Reading targets ... "; cout.flush();
  irtkGreyImage **source = new irtkGreyImage *[numberOfImages];
  for (i = 0; i < numberOfImages; i++) {
      cout << source_name[i] << endl;
      source[i] = new irtkGreyImage(source_name[i]);
  }
  cout << "done" << endl;

  // Create registration filter
  irtkMultipleImageFreeFormRegistration2 *registration = NULL;
  registration = new irtkMultipleImageFreeFormRegistration2;

  // Create initial multi-level free-form deformation
  irtkMultiLevelFreeFormTransformation *mffd = NULL;

  // Default parameters
  padding = MIN_GREY;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-dofout") == 0)) {
      argc--;
      argv++;
      dofout_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-debug") == 0)) {
      argc--;
      argv++;
      ok = true;
      registration->SetDebugFlag(true);
    }

    if ((ok == false) && (strcmp(argv[1], "-x_only") == 0)) {
        argc--;
        argv++;
        registration->SetMode(RegisterX);
        ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-y_only") == 0)) {
        argc--;
        argv++;
        registration->SetMode(RegisterY);
        ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xy_only") == 0)) {
        argc--;
        argv++;
        registration->SetMode(RegisterXY);
        ok = true;
    }

    if ((ok == false) && ((strcmp(argv[1], "-parameter") == 0) || (strcmp(argv[1], "-parin") == 0))) {
      argc--;
      argv++;
      ok = true;
      parin_name = argv[1];
      argc--;
      argv++;
    }
    if ((ok == false) && (strcmp(argv[1], "-parout") == 0)) {
      argc--;
      argv++;
      ok = true;
      parout_name = argv[1];
      argc--;
      argv++;
    }
    if ((ok == false) && (strcmp(argv[1], "-ds") == 0)) {
        argc--;
        argv++;
        spacing = atof(argv[1]);
        argc--;
        argv++;
        ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sp") == 0)) {
        argc--;
        argv++;
        sp = atof(argv[1]);
        argc--;
        argv++;
        ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)) {
      argc--;
      argv++;
      mask_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-mask_dilation") == 0)) {
      argc--;
      argv++;
      mask_dilation = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-weighting") == 0)) {
        argc--;
        argv++;
        int level = atoi(argv[1]);
        if(level > 9) {
            cout << "Number of level < 10" << endl;
            exit(1);
        }
        argc--;
        argv++;
        cout << "Setting weighting for level: " << level << endl; cout.flush();
        for (i = 0; i < numberOfImages; i++) {
            weight[level][i] = atof(argv[1]);
            argc--;
            argv++;
        }
        ok = true;
    }
#ifdef HAS_VTK
    if ((ok == false) && (strcmp(argv[1], "-landmarks") == 0)) {
        argc--;
        argv++;
        cout << "Reading landmark sequence ... "; cout.flush();
        landmarks = new vtkPolyData *[2];
        reader = vtkPolyDataReader::New();
        for (i = 0; i < 2; i++) {
            landmarks[i] = vtkPolyData::New();
            reader->SetFileName(argv[1]);
            reader->Modified();
            reader->Update();
            landmarks[i]->DeepCopy(reader->GetOutput());
            argc--;
            argv++;
        }
        reader->Delete();
        ok = true;
    }
#endif
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Is there a mask to use?
  if (mask_name != NULL) {
    int i;
    irtkGreyImage mask(mask_name);

    if (mask.GetX() != target[0]->GetX() ||
        mask.GetY() != target[0]->GetY() ||
        mask.GetZ() != target[0]->GetZ()) {
      cerr << "Mask given does not match target dimensions." << endl;
      exit(1);
    }

    if (mask_dilation > 0) {
      irtkDilation<irtkGreyPixel> dilation;
    	dilation.SetConnectivity(CONNECTIVITY_26);
      dilation.SetInput(&mask);
      dilation.SetOutput(&mask);
      cout << "Dilating mask ... ";
      cout.flush();
      for (i = 0; i < mask_dilation; i++) dilation.Run();
      cout << "done" << endl;

    }
  }

  if (dofin_name != NULL) {
      irtkTransformation *transform = irtkTransformation::New(dofin_name);
      if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
          mffd = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
      } else {
          if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
              mffd = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
          } else {
              if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
                  mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
              } else {
                  cerr << "Input transformation is not of type rigid or affine " << endl;
                  cerr << "or multi-level free form deformation" << endl;
                  exit(1);
              }
          }
      }
      delete transform;
  } else {
      // Otherwise use identity transformation to start
      mffd = new irtkMultiLevelFreeFormTransformation;
  }

  // Set input and output for the registration filter
  registration->SetInput(target, source, numberOfImages);

#ifdef HAS_VTK
  vtkPolyData* tlandmarks = vtkPolyData::New();
  vtkPolyData* slandmarks = vtkPolyData::New();
  if(landmarks != NULL) {
      tlandmarks->DeepCopy(landmarks[0]);
      slandmarks->DeepCopy(landmarks[1]);
  }
  if(landmarks != NULL) {
      registration->SetLandmarks(tlandmarks,slandmarks);
  }
#endif

  registration->SetOutput(mffd);

  // Set weighting
  registration->SetWeighting(weight);

  // Make an initial Guess for the parameters.
  registration->GuessParameter();
  // Overrride with any the user has set.
  if (parin_name != NULL) {
      registration->irtkMultipleImageRegistration2::Read(parin_name);
  }

  if(sp != 0) registration->SetLambda1(sp);

  // Override parameter settings if necessary
  if (padding != MIN_GREY) {
      registration->SetTargetPadding(padding);
  }
  if (spacing > 0) {
      registration->SetDX(spacing);
      registration->SetDY(spacing);
      registration->SetDZ(spacing);
  }

  // Write parameters if necessary
  if (parout_name != NULL) {
      registration->irtkMultipleImageRegistration2::Write(parout_name);
  }

  // Run registration filter
  registration->Run();

  // Write the final transformation estimate
  if (dofout_name != NULL) {
      mffd->irtkTransformation::Write(dofout_name);
  }
#ifdef HAS_VTK
  tlandmarks->Delete();
  slandmarks->Delete();
#endif
}
