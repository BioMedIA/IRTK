/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkRegistration.h>

// Default filenames
char *source_name = NULL, *target_name = NULL;
char *dofin_name  = NULL, *dofout_name = NULL;
char *parin_name  = NULL, *parout_name = NULL;

void usage()
{
  cerr << "Usage: areg [target] [source] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-parin file>        Read parameter from file" << endl;
  cerr << "<-parout file>       Write parameter to file" << endl;
  cerr << "<-dofin  file>       Read transformation from file" << endl;
  cerr << "<-dofout file>       Write transformation to file" << endl;
  cerr << "<-p3>                Rigid transformation with 3 dofs (for -image)" << endl;
  cerr << "<-p6>                Rigid transformation with 6 dofs (for -image)" << endl;
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
  cerr << "<-image>             Project transformation into image coordinate" << endl;
  cerr << "<-translation_scale> Allow only translation and scale" << endl;
  cerr << "                     before running registration filter." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int ok, padding;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;
  double tox, toy, toz, sox, soy, soz;
  tox = toy = toz = sox = soy = soz = 0.0;
  irtkImageAttributes tatr,satr,atr;
  int centerImages = false, worldImages = false;
  irtkMatrix tmat(4, 4);
  irtkMatrix smat(4, 4);
  irtkMatrix itmat(4, 4);
  irtkMatrix ismat(4, 4);
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
  irtkImageRegistration *registration;

  // Create transformation
  irtkTransformation *transformation = new irtkAffineTransformation;
  // Create registration filter
  if ((target.GetZ() == 1) && (source.GetZ() == 1)) {

    // Registration is 2D
    registration = new irtkImageAffineRegistration2D;

    // Registration is 2D only: Ignore translation along z-axis
    transformation->PutStatus(TZ, _Passive);

    // Registration is 2D only: Ignore rotation round x-axis
    transformation->PutStatus(RX, _Passive);

    // Registration is 2D only: Ignore rotation round y-axis
    transformation->PutStatus(RY, _Passive);

    // Registration is 2D only: Ignore scaling along z-axis
    transformation->PutStatus(SZ, _Passive);

    // Registration is 2D only: Ignore shears involving Z component.
    transformation->PutStatus(SYZ, _Passive);
    transformation->PutStatus(SXZ, _Passive);

  } else {

    // Registration is 3D
    registration = new irtkImageAffineRegistration;

  }

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
	if ((ok == false) && (strcmp(argv[1], "-p3") == 0)) {
      argc--;
      argv++;
	  transformation->PutStatus(RX,  _Passive);
      transformation->PutStatus(RY,  _Passive);
      transformation->PutStatus(RZ,  _Passive);
	  transformation->PutStatus(SX,  _Passive);
      transformation->PutStatus(SY,  _Passive);
      transformation->PutStatus(SZ,  _Passive);
      transformation->PutStatus(SXY, _Passive);
      transformation->PutStatus(SYZ, _Passive);
      transformation->PutStatus(SXZ, _Passive);
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-p6") == 0)) {
		argc--;
		argv++;
		transformation->PutStatus(SX,  _Passive);
		transformation->PutStatus(SY,  _Passive);
		transformation->PutStatus(SZ,  _Passive);
		transformation->PutStatus(SXY, _Passive);
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
	if ((ok == false) && (strcmp(argv[1], "-image") == 0)) {
      argc--;
      argv++;
      worldImages = true;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-translation_scale") == 0)) {
      argc--;
      argv++;
      transformation->PutStatus(RZ,  _Passive);
      transformation->PutStatus(RX,  _Passive);
      transformation->PutStatus(RY,  _Passive);
	  transformation->PutStatus(SXY, _Passive);
      transformation->PutStatus(SYZ, _Passive);
      transformation->PutStatus(SXZ, _Passive);
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
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

  if (worldImages == true) {
    cout << "From world to image ... ";
    // Place the voxel coordinate
	tatr = target.GetImageAttributes();
	satr = source.GetImageAttributes();
	itmat = target.GetImageToWorldMatrix();
	ismat = source.GetWorldToImageMatrix();
	target.PutOrientation(atr._xaxis,atr._yaxis,atr._zaxis);
	//target.PutOrigin(0,0,0);
	target.PutOrigin(double(tatr._x) / 2.0,double(tatr._y) / 2.0,double(tatr._z) / 2.0);
	target.PutPixelSize(atr._dx,atr._dy,atr._dz);
    source.PutOrientation(atr._xaxis,atr._yaxis,atr._zaxis);
	//source.PutOrigin(0,0,0);
	source.PutOrigin(double(satr._x) / 2.0,double(satr._y) / 2.0,double(satr._z) / 2.0);
	source.PutPixelSize(atr._dx,atr._dy,atr._dz);

    irtkAffineTransformation *affineTransf = dynamic_cast<irtkAffineTransformation*> (transformation);

    transfMat = affineTransf->GetMatrix();
    tempMat   = transfMat * itmat;
    tempMat   = ismat * tempMat;

    affineTransf->PutMatrix(tempMat);
    transformation = affineTransf;
    cout << "done" << endl;
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

  // Read parameter if there any, otherwise make an intelligent guess
  if (parin_name != NULL) {
    registration->irtkImageRegistration::Read(parin_name);
  } else {
    registration->GuessParameter();
  }

  if (padding != MIN_GREY) {
    registration->SetTargetPadding(padding);
  }

  // Write parameters if necessary
  if (parout_name != NULL) {
    registration->irtkImageRegistration::Write(parout_name);
  }

  // Run registration filter
  registration->Run();

  // Write the final transformation estimate
  if (dofout_name != NULL) {
    if(centerImages == false && worldImages == false){
		   transformation->Write(dofout_name);
	  }else{
		  irtkAffineTransformation *affineTransf = dynamic_cast<irtkAffineTransformation*> (transformation);
		  transfMat = affineTransf->GetMatrix();
		  tempMat = transfMat;
		  if (centerImages == true) {
			  // Undo the effect of centering the images.	 
			  tmat(0, 3) = -1.0 * tox; 
			  tmat(1, 3) = -1.0 * toy;
			  tmat(2, 3) = -1.0 * toz;
			  smat(0, 3) = sox;
			  smat(1, 3) = soy;
			  smat(2, 3) = soz;

			  tempMat   = tempMat * tmat;
			  tempMat   = smat * tempMat;
		  }
		  if (worldImages == true) {
			  cout << "From image to world ... ";
			  target.PutOrientation(tatr._xaxis,tatr._yaxis,tatr._zaxis);
			  target.PutOrigin(tatr._xorigin,tatr._yorigin,tatr._zorigin);
			  target.PutPixelSize(tatr._dx,tatr._dy,tatr._dz);
			  source.PutOrientation(satr._xaxis,satr._yaxis,satr._zaxis);
			  source.PutOrigin(satr._xorigin,satr._yorigin,satr._zorigin);
			  source.PutPixelSize(satr._dx,satr._dy,satr._dz);
			  itmat = target.GetWorldToImageMatrix();
			  ismat = source.GetImageToWorldMatrix();			 

			  tempMat   = tempMat * itmat;
			  tempMat   = ismat * tempMat;			  
			  cout << "done" << endl;
		  }
		  affineTransf->PutMatrix(tempMat);
		  affineTransf->irtkTransformation::Write(dofout_name);
	  }
  }
}
