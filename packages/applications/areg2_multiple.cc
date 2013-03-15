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

void usage()
{
  cerr << "Usage: areg2_multiple [Number of images] [target1...targetN] [source1...sourceN] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parin file>        Read parameter from file" << endl;
  cerr << "<-parout file>       Write parameter to file" << endl;
  cerr << "<-dofin  file>       Read transformation from file" << endl;
  cerr << "<-dofout file>       Write transformation to file" << endl;
  cerr << "<-p9>                Affine transformation with 9 dofs" << endl;
  cerr << "<-Tp  value>         Padding value in target" << endl;
#ifdef HAS_VTK
  cerr << "<-landmarks t s>     Landmark Regulation target and source" << endl;
#endif
  cerr << "<-center>            Center voxel grids onto image origins" << endl;
  cerr << "                     before running registration." << endl;
  cerr << "<-debug>             Enable debugging information" << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, padding;
  int i, numberOfImages;
  double tox, toy, toz, sox, soy, soz;
  tox = toy = toz = sox = soy = soz = 0.0;
  int centerImages = false;
  irtkMatrix tmat(4, 4);
  irtkMatrix smat(4, 4);
  irtkMatrix tempMat, transfMat;
  tmat.Ident();
  smat.Ident();

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

  // Create transformation
  irtkTransformation *transformation = new irtkAffineTransformation;

  // Create registration
  irtkMultipleImageRegistration2 *registration = new irtkMultipleImageRigidRegistration2;

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
			transformation->PutStatus(TY, _Passive);
			transformation->PutStatus(TZ, _Passive);
			transformation->PutStatus(RX, _Passive);
			transformation->PutStatus(RY, _Passive);
			transformation->PutStatus(RZ, _Passive);
			transformation->PutStatus(SY, _Passive);
			transformation->PutStatus(SZ, _Passive);
			transformation->PutStatus(SXY, _Passive);
			transformation->PutStatus(SYZ, _Passive);
			transformation->PutStatus(SXZ, _Passive);
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-xy_only") == 0)) {
			argc--;
			argv++;
			transformation->PutStatus(TZ, _Passive);
			transformation->PutStatus(RX, _Passive);
			transformation->PutStatus(RY, _Passive);
			transformation->PutStatus(SZ, _Passive);
			transformation->PutStatus(SYZ, _Passive);
			transformation->PutStatus(SXZ, _Passive);
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-p9") == 0)) {
			argc--;
			argv++;
			transformation->PutStatus(SXY, _Passive);
			transformation->PutStatus(SYZ, _Passive);
			transformation->PutStatus(SXZ, _Passive);
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
    if ((ok == false) && (strcmp(argv[1], "-center") == 0)) {
      argc--;
      argv++;
      centerImages = true;
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

  // If there is an initial transformation estimate, read it
  if (dofin_name != NULL) {
    transformation->Read(dofin_name);
  }

  if (centerImages == true) {
    cout << "Centering ... ";
    // Place the voxel centre at the world origin.
    target[0]->GetOrigin(tox, toy, toz);
    source[0]->GetOrigin(sox, soy, soz);
    target[0]->PutOrigin(0.0, 0.0, 0.0);
    source[0]->PutOrigin(0.0, 0.0, 0.0);

    // Update the transformation accordingly.
    tmat(0, 3) = tox;
    tmat(1, 3) = toy;
    tmat(2, 3) = toz;
    smat(0, 3) = -1.0 * sox;
    smat(1, 3) = -1.0 * soy;
    smat(2, 3) = -1.0 * soz;

    irtkAffineTransformation *rigidTransf = dynamic_cast<irtkAffineTransformation*> (transformation);

    transfMat = rigidTransf->GetMatrix();
    tempMat   = transfMat * tmat;
    tempMat   = smat * tempMat;

    rigidTransf->PutMatrix(tempMat);
    transformation = rigidTransf;
    cout << "done" << endl;
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

  registration->SetOutput(transformation);

  // Make an initial Guess for the parameters.
  registration->GuessParameter();
  // Overrride with any the user has set.
  if (parin_name != NULL) {
    registration->irtkMultipleImageRegistration2::Read(parin_name);
  }

  if (padding != MIN_GREY) {
    registration->SetTargetPadding(padding);
  }

  // Write parameters if necessary
  if (parout_name != NULL) {
    registration->irtkMultipleImageRegistration2::Write(parout_name);
  }

  // Run registration filter
  registration->Run();

  // Write the final transformation estimate
  if (dofout_name != NULL) {
    if (centerImages == false) {
      transformation->Write(dofout_name);
    } else {
      // Undo the effect of centering the images.
      irtkAffineTransformation *rigidTransf = dynamic_cast<irtkAffineTransformation*> (transformation);

      tmat(0, 3) = -1.0 * tox;
      tmat(1, 3) = -1.0 * toy;
      tmat(2, 3) = -1.0 * toz;
      smat(0, 3) = sox;
      smat(1, 3) = soy;
      smat(2, 3) = soz;

      transfMat = rigidTransf->GetMatrix();
      tempMat   = transfMat * tmat;
      tempMat   = smat * tempMat;

      rigidTransf->PutMatrix(tempMat);
      rigidTransf->irtkTransformation::Write(dofout_name);
    }
  }
#ifdef HAS_VTK
  tlandmarks->Delete();
  slandmarks->Delete();
#endif
}
