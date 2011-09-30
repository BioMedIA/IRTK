/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#include <irtkCardiac.h>

// Default filenames
char *source_name = NULL, *target_name = NULL;
char *resultinput_name = NULL, *resultout_name = NULL;
char *parin_name  = NULL, *parout_name = NULL;
char **extra_name = NULL, *prefix_name = NULL;

void usage()
{
  cerr << "Usage: spatial_correct [SA target volume]  [SA before correction] [SA after correction] <options> \n" << endl;
  cerr << "where <options> is one or more of the following:\n" << endl;
  cerr << "<-source file>       Read [3D/Detagged image] source from file" << endl;
  cerr << "<-parin file>        Read parameter from file" << endl;
  cerr << "<-parout file>       Write parameter to file" << endl;
  cerr << "<-p3>                Rigid transformation with 3 dofs (for -image)" << endl;
  cerr << "<-p6>                Rigid transformation with 6 dofs (for -image)" << endl;
  cerr << "<-p9>                Affine transformation with 9 dofs" << endl;
  cerr << "<-p12>               Affine transformation with 12 dofs" << endl;
  cerr << "<-Tp  value>         Padding value in target" << endl;
  cerr << "<-debug>             Enable debugging information" << endl;
  cerr << "<-center>            Center voxel grids onto image origins " << endl;
  cerr << "<-translation>       Allow only translation" << endl;
  cerr << "<-extraimage [n] [outputprefix] [input1...inputn]> Allow extra images transformed" << endl;
  cerr << "                     before running registration filter." << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, t, ok, padding, number_extraImage;
  int target_x1, target_y1, target_z1, target_x2, target_y2, target_z2;
  int source_x1, source_y1, source_z1, source_x2, source_y2, source_z2;
  double tox, toy, toz, sox, soy, soz;
  tox = toy = toz = sox = soy = soz = 0.0;
  irtkImageAttributes tatr,satr,atr;
  int centerImages = false, worldImages = true;
  irtkMatrix tmat(4, 4);
  irtkMatrix smat(4, 4);
  irtkMatrix itmat(4, 4);
  irtkMatrix ismat(4, 4);
  irtkMatrix tempMat, transfMat;
  tmat.Ident();
  smat.Ident();

  // Check command line
  if (argc < 4) {
    usage();
  }

  // Parse source and target images
  target_name = argv[1];
  argc--;
  argv++;
  resultinput_name = argv[1];
  argc--;
  argv++;
  resultout_name = argv[1];
  argc--;
  argv++;

  // Read target image
  cout << "Reading target ... "; cout.flush();
  irtkGreyImage target(target_name);
  cout << "done" << endl;
  // Read input image
  cout << "Reading SA ... "; cout.flush();
  irtkGreyImage input(resultinput_name);
  cout << "done" << endl;
  irtkGreyImage output(input.GetImageAttributes());

  // Create registration filter
  irtkCardiacSpatialCorrection *registration;

  // Create transformation
  irtkTransformation **transformation = new irtkTransformation*[target.GetZ()];

  for(k = 0; k < target.GetZ(); k++){
	  transformation[k] = new irtkAffineTransformation;
  }

  // Registration is 3D
  registration = new irtkCardiacSpatialCorrection;

  // Fix ROI
  target_x1 = 0;
  target_y1 = 0;
  target_z1 = 0;
  target_x2 = target.GetX();
  target_y2 = target.GetY();
  target_z2 = target.GetZ();

  // Extra image number
  number_extraImage = 0;

  // Default parameters
  padding = MIN_GREY;

  // Parse remaining parameters
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-Tp") == 0)) {
      argc--;
      argv++;
      padding = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
	if ((ok == false) && ((strcmp(argv[1], "-source") == 0))) {
		argc--;
		argv++;
		ok = true;
		source_name = argv[1];
		argc--;
		argv++;
	}
    if ((ok == false) && (strcmp(argv[1], "-x_only") == 0)) {
      argc--;
      argv++;
	  for(k = 0; k < target.GetZ(); k++){
		  transformation[k]->PutStatus(TY,  _Passive);
		  transformation[k]->PutStatus(TZ,  _Passive);
		  transformation[k]->PutStatus(RX,  _Passive);
		  transformation[k]->PutStatus(RY,  _Passive);
		  transformation[k]->PutStatus(RZ,  _Passive);
		  transformation[k]->PutStatus(SY,  _Passive);
		  transformation[k]->PutStatus(SZ,  _Passive);
		  transformation[k]->PutStatus(SXY, _Passive);
		  transformation[k]->PutStatus(SYZ, _Passive);
		  transformation[k]->PutStatus(SXZ, _Passive);
	  }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-xy_only") == 0)) {
      argc--;
      argv++;
	  for(k = 0; k < target.GetZ(); k++){
		  transformation[k]->PutStatus(TZ,  _Passive);
		  transformation[k]->PutStatus(RX,  _Passive);
		  transformation[k]->PutStatus(RY,  _Passive);
		  transformation[k]->PutStatus(SZ,  _Passive);
		  transformation[k]->PutStatus(SYZ, _Passive);
		  transformation[k]->PutStatus(SXZ, _Passive);
	  }
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-p3") == 0)) {
      argc--;
      argv++;
	  for(k = 0; k < target.GetZ(); k++){
		  transformation[k]->PutStatus(RX,  _Passive);
		  transformation[k]->PutStatus(RY,  _Passive);
		  transformation[k]->PutStatus(RZ,  _Passive);
		  transformation[k]->PutStatus(SX,  _Passive);
		  transformation[k]->PutStatus(SY,  _Passive);
		  transformation[k]->PutStatus(SZ,  _Passive);
		  transformation[k]->PutStatus(SXY, _Passive);
		  transformation[k]->PutStatus(SYZ, _Passive);
		  transformation[k]->PutStatus(SXZ, _Passive);
	  }
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-p6") == 0)) {
		argc--;
		argv++;
		for(k = 0; k < target.GetZ(); k++){
		transformation[k]->PutStatus(SX,  _Passive);
		transformation[k]->PutStatus(SY,  _Passive);
		transformation[k]->PutStatus(SZ,  _Passive);
		transformation[k]->PutStatus(SXY, _Passive);
		transformation[k]->PutStatus(SYZ, _Passive);
		transformation[k]->PutStatus(SXZ, _Passive);
		}
		ok = true;
	}
    if ((ok == false) && (strcmp(argv[1], "-p9") == 0)) {
      argc--;
      argv++;
	  for(k = 0; k < target.GetZ(); k++){
      transformation[k]->PutStatus(SXY, _Passive);
      transformation[k]->PutStatus(SYZ, _Passive);
      transformation[k]->PutStatus(SXZ, _Passive);
	  }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-p12") == 0)) {
      argc--;
      argv++;
      ok = true;
    }
	if ((ok == false) && (strcmp(argv[1], "-extraimage") == 0)) {
		argc--;
		argv++;
		number_extraImage = atoi(argv[1]);
		argc--;
		argv++;
		prefix_name = argv[1];
		argc--;
		argv++;
		extra_name = new char*[number_extraImage];
		// Read atlas for each tissue
		for (i = 0; i < number_extraImage; i++) {
			extra_name[i] = argv[1];
			cerr << "Image " << i <<" = " << argv[1] <<endl;
			argc--;
			argv++;
		}
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
	if ((ok == false) && (strcmp(argv[1], "-translation") == 0)) {
      argc--;
      argv++;
	  for(k = 0; k < target.GetZ(); k++){
		  transformation[k]->PutStatus(RZ,  _Passive);
		  transformation[k]->PutStatus(RX,  _Passive);
		  transformation[k]->PutStatus(RY,  _Passive);
		  transformation[k]->PutStatus(SX,  _Passive);
		  transformation[k]->PutStatus(SY,  _Passive);
		  transformation[k]->PutStatus(SZ,  _Passive);
		  transformation[k]->PutStatus(SXY, _Passive);
		  transformation[k]->PutStatus(SYZ, _Passive);
		  transformation[k]->PutStatus(SXZ, _Passive);
	  }
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkGreyImage source;

  // Read source image
  if(source_name == NULL){
	  cout << "Creating source ... "; cout.flush();
	  source.Initialize(target.GetImageAttributes());
	  cout << "done" << endl;
  }
  else{
	  cout << "Reading source ... "; cout.flush();
	  source.Read(source_name);
	  cout << "done" << endl;
  }

  source_x1 = 0;
  source_y1 = 0;
  source_z1 = 0;
  source_x2 = source.GetX();
  source_y2 = source.GetY();
  source_z2 = source.GetZ();

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

  if (worldImages == true) {
    cout << "From world to image ... ";
    // Place the voxel coordinate
	tatr = target.GetImageAttributes();
	satr = source.GetImageAttributes();
	itmat = target.GetImageToWorldMatrix();
	ismat = source.GetWorldToImageMatrix();
	target.PutOrientation(atr._xaxis,atr._yaxis,atr._zaxis);
	target.PutOrigin(double(tatr._x) / 2.0,double(tatr._y) / 2.0,double(tatr._z) / 2.0);
	target.PutPixelSize(atr._dx,atr._dy,atr._dz);
    source.PutOrientation(atr._xaxis,atr._yaxis,atr._zaxis);
	source.PutOrigin(double(satr._x) / 2.0,double(satr._y) / 2.0,double(satr._z) / 2.0);
	source.PutPixelSize(atr._dx,atr._dy,atr._dz);

	for(k = 0; k < target.GetZ(); k++){
		irtkAffineTransformation *affineTransf = dynamic_cast<irtkAffineTransformation*> (transformation[k]);

		transfMat = affineTransf->GetMatrix();
		tempMat   = transfMat * itmat;
		tempMat   = ismat * tempMat;

		affineTransf->PutMatrix(tempMat);
		transformation[k] = affineTransf;
	}
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

	for(k = 0; k < target.GetZ(); k++){
		irtkAffineTransformation *affineTransf = dynamic_cast<irtkAffineTransformation*> (transformation[k]);

		transfMat = affineTransf->GetMatrix();
		tempMat   = transfMat * tmat;
		tempMat   = smat * tempMat;

		affineTransf->PutMatrix(tempMat);
		transformation[k] = affineTransf;
	}
    cout << "done" << endl;
  }

  // Set input and output for the registration filter
  registration->SetInput(&target, &source);
  registration->SetOutput(transformation);

  // Make an initial Guess for the parameters.
  registration->GuessParameter();
  // Overrride with any the user has set.
  if (parin_name != NULL) {
    registration->Read(parin_name);
  }

  if (padding != MIN_GREY) {
    registration->SetTargetPadding(padding);
  }

  // Write parameters if necessary
  if (parout_name != NULL) {
    registration->Write(parout_name);
  }

  // Run registration filter
  registration->Run();

  // Correct the final transformation estimate
  if(centerImages == false && worldImages == false){

  }else{
	  for(k = 0; k < target.GetZ(); k++){
		  irtkAffineTransformation *affineTransf = dynamic_cast<irtkAffineTransformation*> (transformation[k]);
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
		  transformation[k] = affineTransf;
	  }
	}
  //write the final output

  for(k = 0; k < input.GetZ(); k++){
	  //invert transformation
	  ((irtkAffineTransformation*)transformation[k])->Invert();
	  for(j = 0; j < input.GetY(); j++){
		  for(i = 0; i < input.GetX(); i++){
			  //transformation
			  double x,y,z;
			  x = i; y = j; z = k;
			  output.ImageToWorld(x,y,z);
			  transformation[k]->Transform(x,y,z);
			  input.WorldToImage(x,y,z);
			  for(t = 0; t < input.GetT(); t++){
				  if(round(x)>=0 && round(x)<input.GetX()
					  && round(y)>=0 && round(y)<input.GetY()
					  && round(z)>=0 && round(z)<input.GetZ()){
						  output.PutAsDouble(i,j,k,t,
							  input.GetAsDouble(round(x),round(y),round(z),t));
				  }else{
					  output.PutAsDouble(i,j,k,t,0);
				  }
			  }
		  }
	  }
  }
  output.Write(resultout_name);

  if(number_extraImage != 0){
	  for(int n = 0; n < number_extraImage; n++){
		  irtkGreyImage extrainput(extra_name[n]);
		  irtkGreyImage extraoutput(extrainput.GetImageAttributes());
		  for(k = 0; k < extrainput.GetZ(); k++){
			  //invert transformation
			  for(j = 0; j < extrainput.GetY(); j++){
				  for(i = 0; i < extrainput.GetX(); i++){
					  //transformation
					  double x,y,z;
					  x = i; y = j; z = k;
					  extraoutput.ImageToWorld(x,y,z);
					  transformation[k]->Transform(x,y,z);
					  extrainput.WorldToImage(x,y,z);
					  for(t = 0; t < extrainput.GetT(); t++){
						  if(round(x)>=0 && round(x)<extrainput.GetX()
							  && round(y)>=0 && round(y)<extrainput.GetY()
							  && round(z)>=0 && round(z)<extrainput.GetZ()){
								  extraoutput.PutAsDouble(i,j,k,t,
									  extrainput.GetAsDouble(round(x),round(y),round(z),t));
						  }else{
							  extraoutput.PutAsDouble(i,j,k,t,0);
						  }
					  }
				  }
			  }
		  }
		  char buffer[255];
		  sprintf(buffer,"%s%s",prefix_name,extra_name[n]);
		  extraoutput.Write(buffer);
	  }
  }

  for(k = 0; k < target.GetZ(); k++){
	  delete transformation[k];
  }

  delete []transformation;
}
