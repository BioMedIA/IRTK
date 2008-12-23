/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkDemonsRegistration.h>

// Default filenames
char *source_name = NULL, *target_name = NULL;
char *dofin_name  = NULL, *dofout_name = NULL;
char *parameter_name = NULL;

int main(int argc, char **argv)
{
  int ok;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;
  double dx, dy, dz;

  // Create registration filter
  irtkDemonsRegistration *registration = new irtkDemonsRegistration;

  // Create initial multi-level free-form deformation
  irtkMultiLevelFreeFormTransformation *mffd = NULL;

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  source_name = argv[1];
  argc--;
  argv++;

  // Read target image
  cout << "Reading target ... "; cout.flush();
  irtkGreyImage target(target_name);
  cout << "done" << endl;
  // Read source image
  cout << "Reading source ... "; cout.flush();
  irtkGreyImage source(source_name);
  cout << "done" << endl;

  // Fix ROI
  target_x1 = 0;
  target_y1 = 0;
  target_z1 = 0;
  target_x2 = target.GetX();
  target_y2 = target.GetY();
  target_z2 = target.GetZ();
  source_x1 = 0;
  source_y1 = 0;
  source_z1 = 0;
  source_x2 = source.GetX();
  source_y2 = source.GetY();
  source_z2 = source.GetZ();

  // Parse remaining parameters
  dx = 1;
  dy = 1;
  dz = 1;
  while (argc > 1) {
    ok = False;
    if ((ok == False) && (strcmp(argv[1], "-dofin") == 0)) {
      argc--;
      argv++;
      dofin_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-dofout") == 0)) {
      argc--;
      argv++;
      dofout_name = argv[1];
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      target_x1 = source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      target_x2 = source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      target_y1 = source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      target_y2 = source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      target_z1 = source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      target_z2 = source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tx1") == 0)) {
      argc--;
      argv++;
      target_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tx2") == 0)) {
      argc--;
      argv++;
      target_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ty1") == 0)) {
      argc--;
      argv++;
      target_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Ty2") == 0)) {
      argc--;
      argv++;
      target_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tz1") == 0)) {
      argc--;
      argv++;
      target_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tz2") == 0)) {
      argc--;
      argv++;
      target_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sx1") == 0)) {
      argc--;
      argv++;
      source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sx2") == 0)) {
      argc--;
      argv++;
      source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sy1") == 0)) {
      argc--;
      argv++;
      source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sy2") == 0)) {
      argc--;
      argv++;
      source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sz1") == 0)) {
      argc--;
      argv++;
      source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Sz2") == 0)) {
      argc--;
      argv++;
      source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      registration->SetTargetPadding(atoi(argv[1]));
      argc--;
      argv++;
      ok = True;
    }
    if ((ok == False) && (strcmp(argv[1], "-parameter") == 0)) {
      argc--;
      argv++;
      ok = True;
      parameter_name = argv[1];
      argc--;
      argv++;
    }
    if (ok == False) {
      cerr << "Can not parse argument " << argv[1] << endl;
    }
  }

  // If there is an region of interest, use it
  if ((target_x1 != 0) || (target_x2 != target.GetX()) ||
      (target_y1 != 0) || (target_y2 != target.GetY()) ||
      (target_z1 != 0) || (target_z2 != target.GetZ())) {
    target = target.GetRegion(target_x1, target_y1, target_z1,
                              target_x2, target_y2, target_z2);
  }

  // If there is an region of interest for the source image, use it
  if ((source_x1 != 0) || (source_x2 != source.GetX()) ||
      (source_y1 != 0) || (source_y2 != source.GetY()) ||
      (source_z1 != 0) || (source_z2 != source.GetZ())) {
    source = source.GetRegion(source_x1, source_y1, source_z1,
                              source_x2, source_y2, source_z2);
  }

  if (dofin_name != NULL) {
    if (irtkRigidTransformation::CheckHeader(dofin_name) == True) {
      mffd = new irtkMultiLevelFreeFormTransformation;
      mffd->irtkTransformation::Read(dofin_name);
    } else {
      if (irtkAffineTransformation::CheckHeader(dofin_name) == True) {
        mffd = new irtkMultiLevelFreeFormTransformation;
        mffd->irtkTransformation::Read(dofin_name);
      } else {
        if (irtkMultiLevelFreeFormTransformation::CheckHeader(dofin_name) == True) {
          mffd = new irtkMultiLevelFreeFormTransformation;
          mffd->irtkTransformation::Read(dofin_name);
        } else {
          cerr << "Input transformation is not of type rigid, affine " << endl;
          cerr << "or multi-level free form deformation" << endl;
          exit(1);
        }
      }
    }
  } else {
    // Otherwise use identity transformation to start
    mffd = new irtkMultiLevelFreeFormTransformation;
  }

  // and add ffd to multi-level ffd
  irtkLinearFreeFormTransformation *ffd = new irtkLinearFreeFormTransformation(target, dx, dy, dz);

  // Add ffd
  mffd->PushLocalTransformation(ffd);

  // Set input and output for the registration filter
  registration->SetInput(&target, &source);
  registration->SetOutput(mffd);

  // Read parameter if there any, otherwise make an intelligent guess
  if (parameter_name != NULL) {
    registration->Read(parameter_name);
  }

  // Run registration filter
  registration->Run();

  // Write the final transformation estimate
  if (dofout_name != NULL) {
    mffd->irtkTransformation::Write(dofout_name);
  }
}
