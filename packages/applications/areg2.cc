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
char *source_name = NULL, *target_name = NULL;
char *dofin_name  = NULL, *dofout_name = NULL;
char *parin_name  = NULL, *parout_name = NULL;
char *mask_name = NULL;

void usage()
{
  cerr << "Usage: areg2 [target] [source] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parin file>        Read parameter from file" << endl;
  cerr << "<-parout file>       Write parameter to file" << endl;
  cerr << "<-dofin  file>       Read transformation from file" << endl;
  cerr << "<-dofout file>       Write transformation to file" << endl;
  cerr << "<-p9>                Affine transformation with 9 dofs" << endl;
  cerr << "<-p12>               Affine transformation with 12 dofs" << endl;
  cerr << "<-Rx1 value>         Region of interest in both images" << endl;
  cerr << "<-Ry1 value>         Region of interest in both images" << endl;
  cerr << "<-Rz1 value>         Region of interest in both images" << endl;
  cerr << "<-Rx2 value>         Region of interest in both images" << endl;
  cerr << "<-Ry2 value>         Region of interest in both images" << endl;
  cerr << "<-Rz2 value>         Region of interest in both images" << endl;
  cerr << "<-Tx1 value>         Region of interest in target image" << endl;
  cerr << "<-Ty1 value>         Region of interest in target image" << endl;
  cerr << "<-Tz1 value>         Region of interest in target image" << endl;
  cerr << "<-Tx2 value>         Region of interest in target image" << endl;
  cerr << "<-Ty2 value>         Region of interest in target image" << endl;
  cerr << "<-Tz2 value>         Region of interest in target image" << endl;
  cerr << "<-Sx1 value>         Region of interest in source image" << endl;
  cerr << "<-Sy1 value>         Region of interest in source image" << endl;
  cerr << "<-Sz1 value>         Region of interest in source image" << endl;
  cerr << "<-Sx2 value>         Region of interest in source image" << endl;
  cerr << "<-Sy2 value>         Region of interest in source image" << endl;
  cerr << "<-Sz2 value>         Region of interest in source image" << endl;
  cerr << "<-Tp  value>         Padding value in target" << endl;
  cerr << "<-debug>             Enable debugging information" << endl;
  cerr << "<-center>            Center voxel grids onto image origins " << endl;
  cerr << "                     before running registration filter." << endl;
  cerr << "<-mask file>         Use a mask to define the ROI. The mask" << endl;
  cerr << "                     must have the same dimensions as the target." << endl;
  cerr << "                     Voxels in the mask with zero or less are " << endl;
  cerr << "                     padded in the target." << endl;
  cerr << "<-mask_dilation n>   Dilate mask n times before using it" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int ok, padding;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;
  double tox, toy, toz, sox, soy, soz;
  tox = toy = toz = sox = soy = soz = 0.0;
  int centerImages = false, mask_dilation = 0;
  irtkMatrix tmat(4, 4);
  irtkMatrix smat(4, 4);
  irtkMatrix tempMat, transfMat;
  tmat.Ident();
  smat.Ident();

  // Check command line
  if (argc < 3) {
    usage();
  }

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

  // Create registration filter
  irtkImageRegistration2 *registration;

  // Create transformation
  irtkTransformation *transformation = new irtkAffineTransformation;

  // Registration is 3D
  registration = new irtkImageAffineRegistration2;

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
    if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
      argc--;
      argv++;
      target_x1 = source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
      argc--;
      argv++;
      target_x2 = source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
      argc--;
      argv++;
      target_y1 = source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
      argc--;
      argv++;
      target_y2 = source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
      argc--;
      argv++;
      target_z1 = source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
      argc--;
      argv++;
      target_z2 = source_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tx1") == 0)) {
      argc--;
      argv++;
      target_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tx2") == 0)) {
      argc--;
      argv++;
      target_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty1") == 0)) {
      argc--;
      argv++;
      target_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Ty2") == 0)) {
      argc--;
      argv++;
      target_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz1") == 0)) {
      argc--;
      argv++;
      target_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Tz2") == 0)) {
      argc--;
      argv++;
      target_z2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sx1") == 0)) {
      argc--;
      argv++;
      source_x1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sx2") == 0)) {
      argc--;
      argv++;
      source_x2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sy1") == 0)) {
      argc--;
      argv++;
      source_y1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sy2") == 0)) {
      argc--;
      argv++;
      source_y2 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sz1") == 0)) {
      argc--;
      argv++;
      source_z1 = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Sz2") == 0)) {
      argc--;
      argv++;
      source_z2 = atoi(argv[1]);
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
    if ((ok == false) && (strcmp(argv[1], "-x_only") == 0)) {
      argc--;
      argv++;
      transformation->PutStatus(TY,  _Passive);
      transformation->PutStatus(TZ,  _Passive);
      transformation->PutStatus(RX,  _Passive);
      transformation->PutStatus(RY,  _Passive);
      transformation->PutStatus(RZ,  _Passive);
      transformation->PutStatus(SY,  _Passive);
      transformation->PutStatus(SZ,  _Passive);
      transformation->PutStatus(SXY, _Passive);
      transformation->PutStatus(SYZ, _Passive);
      transformation->PutStatus(SXZ, _Passive);
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xy_only") == 0)) {
      argc--;
      argv++;
      transformation->PutStatus(TZ,  _Passive);
      transformation->PutStatus(RX,  _Passive);
      transformation->PutStatus(RY,  _Passive);
      transformation->PutStatus(SZ,  _Passive);
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
    if ((ok == false) && (strcmp(argv[1], "-p12") == 0)) {
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
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Is there a mask to use?
  if (mask_name != NULL) {
    int voxels, i;
    irtkGreyPixel *ptr2target, *ptr2mask;
    irtkGreyImage mask(mask_name);

    if (mask.GetX() != target.GetX() ||
        mask.GetY() != target.GetY() ||
        mask.GetZ() != target.GetZ()) {
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

    voxels     = target.GetNumberOfVoxels();
    ptr2target = target.GetPointerToVoxels();
    ptr2mask   = mask.GetPointerToVoxels();

    for (i = 0; i < voxels; ++i) {
      if (*ptr2mask <= 0) {
        *ptr2target = padding;
      }
      ++ptr2mask;
      ++ptr2target;
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

  // If there is an initial transformation estimate, read it
  if (dofin_name != NULL) {
    transformation->Read(dofin_name);
  }

  if (centerImages == true) {
    cout << "Centering ... ";
    // Place the voxel centre at the world origin.
    target.GetOrigin(tox, toy, toz);
    source.GetOrigin(sox, soy, soz);
    target.PutOrigin(0.0, 0.0, 0.0);
    source.PutOrigin(0.0, 0.0, 0.0);

    // Update the transformation accordingly.
    tmat(0, 3) = tox;
    tmat(1, 3) = toy;
    tmat(2, 3) = toz;
    smat(0, 3) = -1.0 * sox;
    smat(1, 3) = -1.0 * soy;
    smat(2, 3) = -1.0 * soz;

    irtkAffineTransformation *affineTransf = dynamic_cast<irtkAffineTransformation*> (transformation);

    transfMat = affineTransf->GetMatrix();
    tempMat   = transfMat * tmat;
    tempMat   = smat * tempMat;

    affineTransf->PutMatrix(tempMat);
    transformation = affineTransf;
    cout << "done" << endl;
  }

  // Set input and output for the registration filter
  registration->SetInput(&target, &source);
  registration->SetOutput(transformation);

  // Make an initial Guess for the parameters.
  registration->GuessParameter();
  // Overrride with any the user has set.
  if (parin_name != NULL) {
    registration->irtkImageRegistration2::Read(parin_name);
  }

  if (padding != MIN_GREY) {
    registration->SetTargetPadding(padding);
  }

  // Write parameters if necessary
  if (parout_name != NULL) {
    registration->irtkImageRegistration2::Write(parout_name);
  }

  // Run registration filter
  registration->Run();

  // Write the final transformation estimate
  if (dofout_name != NULL) {
    if (centerImages == false) {
      transformation->Write(dofout_name);
    } else {
      // Undo the effect of centering the images.
      irtkAffineTransformation *affineTransf = dynamic_cast<irtkAffineTransformation*> (transformation);

      tmat(0, 3) = -1.0 * tox;
      tmat(1, 3) = -1.0 * toy;
      tmat(2, 3) = -1.0 * toz;
      smat(0, 3) = sox;
      smat(1, 3) = soy;
      smat(2, 3) = soz;

      transfMat = affineTransf->GetMatrix();
      tempMat   = transfMat * tmat;
      tempMat   = smat * tempMat;

      affineTransf->PutMatrix(tempMat);
      affineTransf->irtkTransformation::Write(dofout_name);
    }
  }
}
