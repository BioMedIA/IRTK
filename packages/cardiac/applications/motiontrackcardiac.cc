#include <irtkCardiac.h>

#include <irtkGaussianBlurring.h>

#include <irtkGaussianBlurring4D.h>

char *dofout_name = NULL, *parin_name  = NULL, *parin_name2  = NULL, *parout_name = NULL, **untagfilenames = NULL,**tagfilenames = NULL, *thresholdname = NULL;

void numberofframeerror(){
	cerr << "the images' number of frames do not equal" << endl;
	exit(1);
}

void usage()
{
	cerr << "Usage: motiontrackcardiac numberOfUntaggedImages numberOfTaggedImages [untaggedimage sequence] [taggedimage sequence] <options>" << endl;
	cerr << "Registration using all image sequences at the same time and combine similarity measure based on spatially adaptive weight" << endl;
    cerr << "Recommended parameter setting is Lregulation = 0.02 Lambda2 = 0.8 and Lambda1 = 0.0001, -adaptive 0.9" << endl;
    cerr << "result is highly sensitive to parameters due to the complex nature of the algorithm" << endl;
	cerr << "where <options> is one or more of the following:" << endl;
	cerr << "<-threshold file>    Read segmentation from file (must option),myocardium 3" << endl;
	cerr << "<-parin file>        Read parameter from file" << endl;
	cerr << "<-parin2 file>       Read parameter from file for temporal misalignment cases" << endl;
	cerr << "<-parout file>       Write parameter to file" << endl;
	cerr << "<-dofout folder>     Write transformations to folder" << endl;
	cerr << "<-ref file>          Reference time frame (default = first frame of uimage sequence)" << endl;
	cerr << "<-Rx1 value>         Region of interest in images" << endl;
	cerr << "<-Ry1 value>         Region of interest in images" << endl;
	cerr << "<-Rz1 value>         Region of interest in images" << endl;
	cerr << "<-Rt1 value>         Region of interest in images" << endl;
	cerr << "<-Rx2 value>         Region of interest in images" << endl;
	cerr << "<-Ry2 value>         Region of interest in images" << endl;
	cerr << "<-Rz2 value>         Region of interest in images" << endl;
	cerr << "<-Rt2 value>         Region of interest in images" << endl;
	cerr << "<-Tp  value>         Padding value" << endl;
	cerr << "<-landmarks name>    Landmark Regulation input name is prefix" << endl;
    cerr << "<-adaptive weight>   Adapt regulation weight using the weight" <<endl;
	exit(1);
}


int main(int argc, char **argv)
{
	int l, i, numberOfUntaggedImages, numberOfTaggedImages, t, x, y, z, 
		ux1, uy1, uz1, ut1, ux2, uy2, uz2, ut2, ok, debug
		,tx1, ty1, tz1, tt1, tx2, ty2, tz2, tt2, times;
	double spacing, sigma, adaptive, weight;
	irtkGreyPixel padding;
	irtkMultiLevelFreeFormTransformation *mffd;
	irtkGreyImage *threshold = NULL;
	vtkPolyData **landmarks;
    vtkPolyDataReader *reader;
	// Check command line
	if (argc < 6) {
		usage();
	}

	// Get uimage names for sequence
	landmarks = NULL;
    reader = NULL;
	numberOfUntaggedImages = 0;
	numberOfTaggedImages = 0;
	numberOfUntaggedImages = atoi(argv[1]);
	argv++;
	argc--;
	numberOfTaggedImages = atoi(argv[1]);
	argv++;
	argc--;
	untagfilenames = argv;
	untagfilenames++;
	for(i=0;i<numberOfUntaggedImages;i++){
		if ((argc > 1) && (argv[1][0] != '-' )) {
			argv++;
			argc--;
		}else{
			usage();
		}
	}
	tagfilenames = argv;
	tagfilenames++;
	for(i=0;i<numberOfTaggedImages;i++){
		if ((argc > 1) && (argv[1][0] != '-' )) {
			argv++;
			argc--;
		}else{
			usage();
		}
	}

	// Read uimage sequence
	cout << "Reading untagged uimage sequences ... "; cout.flush();
	irtkGreyImage **uimage = new irtkGreyImage *[numberOfUntaggedImages];
	for (i = 0; i < numberOfUntaggedImages; i++) {
		cout << untagfilenames[i] << endl;
		uimage[i] = new irtkGreyImage(untagfilenames[i]);
	}
	cout << "Reading tagged uimage sequences ... "; cout.flush();
	irtkGreyImage **timage = new irtkGreyImage *[numberOfTaggedImages];
	for (i = 0; i < numberOfTaggedImages; i++) {
		cout << tagfilenames[i] << endl;
		timage[i] = new irtkGreyImage(tagfilenames[i]);
	}

	// Check if there is at least one uimage
	if (numberOfTaggedImages == 0 || numberOfUntaggedImages == 0) {
		usage();
	}

	// Fix ROI
	ux1 = 0;
	uy1 = 0;
	uz1 = 0;
	ut1 = 0;
	ux2 = uimage[0]->GetX();
	uy2 = uimage[0]->GetY();
	uz2 = uimage[0]->GetZ();
	ut2 = uimage[0]->GetT();

	tx1 = 0;
	ty1 = 0;
	tz1 = 0;
	tt1 = 0;
	tx2 = timage[0]->GetX();
	ty2 = timage[0]->GetY();
	tz2 = timage[0]->GetZ();
	tt2 = timage[0]->GetT();

	// Default parameters
	padding   = MIN_GREY;
	spacing   = 0;
	sigma     = 0;
    adaptive  = 0;
	debug     = false;

	// Parse remaining parameters
	while (argc > 1) {
		ok = false;
		if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
			argc--;
			argv++;
			ux1 = atoi(argv[1]);
			tx1 = ux1;
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
			argc--;
			argv++;
			ux2 = atoi(argv[1]);
			tx2 = ux2;
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
			argc--;
			argv++;
			uy1 = atoi(argv[1]);
			ty1 = uy1;
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
			argc--;
			argv++;
			uy2 = atoi(argv[1]);
			ty2 = uy2;
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
			argc--;
			argv++;
			uz1 = atoi(argv[1]);
			tz1 = uz1;
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
			argc--;
			argv++;
			uz2 = atoi(argv[1]);
			tz2 = uz2;
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rt1") == 0)) {
			argc--;
			argv++;
			ut1 = atoi(argv[1]);
			tt1 = ut1;
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-Rt2") == 0)) {
			argc--;
			argv++;
			ut2 = atoi(argv[1]);
			tt2 = ut2;
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
		if ((ok == false) && (strcmp(argv[1], "-ds") == 0)) {
			argc--;
			argv++;
			spacing = atof(argv[1]);
			argc--;
			argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-parin") == 0)) {
			argc--;
			argv++;
			ok = true;
			parin_name = argv[1];
			argc--;
			argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-parin2") == 0)) {
			argc--;
			argv++;
			ok = true;
			parin_name2 = argv[1];
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
		if ((ok == false) && (strcmp(argv[1], "-blur") == 0)) {
			argc--;
			argv++;
			ok = true;
			sigma = atof(argv[1]);
			argc--;
			argv++;
		}
        if ((ok == false) && (strcmp(argv[1], "-adaptive") == 0)) {
            argc--;
            argv++;
            ok = true;
            adaptive = atof(argv[1]);
            argc--;
            argv++;
        }
		if ((ok == false) && (strcmp(argv[1], "-threshold") == 0)) {
			argc--;
			argv++;
			ok = true;
			thresholdname = argv[1];
			argc--;
			argv++;
		}
		if ((ok == false) && (strcmp(argv[1], "-debug") == 0)) {
			argc--;
			argv++;
			ok = true;
			debug = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-landmarks") == 0)) {
            argc--;
            argv++;
            cout << "Reading landmark sequence ... "; cout.flush();
            landmarks = new vtkPolyData *[uimage[0]->GetT()];
            reader = vtkPolyDataReader::New();
            for (i = 0; i < uimage[0]->GetT(); i++) {
                char buffer[255];
                sprintf(buffer, "%s%.2d.vtk", argv[1],i);
                cout << buffer << endl;
                landmarks[i] = vtkPolyData::New();            
                reader->SetFileName(buffer);
                reader->Modified();
                reader->Update();
                landmarks[i]->DeepCopy(reader->GetOutput());
            }
            reader->Delete();
            argc--;
            argv++;
            ok = true;
		}
		if (ok == false) {
			cerr << "Can not parse argument " << argv[1] << endl;
			usage();
		}
	}

	if(thresholdname!=NULL)
		threshold = new irtkGreyImage(thresholdname);
	else{
		cerr<<"please generate threshold using graphcut or ems or whatever"<<endl;
        exit(1);
    }

    if(uimage[0]->GetT() == 1){
        cerr << "image only has one frame, can't do motion track" << endl;
        exit(1);
    }

	// If there is an region of interest, use it
	if ((ux1 != 0) || (ux2 != uimage[0]->GetX()) ||
		(uy1 != 0) || (uy2 != uimage[0]->GetY()) ||
		(uz1 != 0) || (uz2 != uimage[0]->GetZ()) ||
		(ut1 != 0) || (ut2 != uimage[0]->GetT())) {
			for (i = 0; i < numberOfUntaggedImages; i++) {
				*uimage[i] = uimage[i]->GetRegion(ux1, uy1, uz1, ut1, ux2, uy2, uz2, ut2);
			}
	}

	// If there is an region of interest, use it
	if ((tx1 != 0) || (tx2 != timage[0]->GetX()) ||
		(ty1 != 0) || (ty2 != timage[0]->GetY()) ||
		(tz1 != 0) || (tz2 != timage[0]->GetZ()) ||
		(tt1 != 0) || (tt2 != timage[0]->GetT())) {
			for (i = 0; i < numberOfTaggedImages; i++) {
				*timage[i] = timage[i]->GetRegion(tx1, ty1, tz1, tt1, tx2, ty2, tz2, tt2);
			}
	}

	// If sigma is larger than 0, blur images using 4D blurring
	if (sigma > 0) {
		cout << "Blurring uimage sequences ... "; cout.flush();
		for (i = 0; i < numberOfUntaggedImages; i++) {
			irtkGaussianBlurring4D<irtkGreyPixel> gaussianBlurring4D(sigma);
			gaussianBlurring4D.SetInput (uimage[i]);
			gaussianBlurring4D.SetOutput(uimage[i]);
			gaussianBlurring4D.Run();
		}
		cout << "done" << endl;
	}
	if (sigma > 0) {
		cout << "Blurring timage sequences ... "; cout.flush();
		for (i = 0; i < numberOfTaggedImages; i++) {
			irtkGaussianBlurring4D<irtkGreyPixel> gaussianBlurring4D(sigma);
			gaussianBlurring4D.SetInput (timage[i]);
			gaussianBlurring4D.SetOutput(timage[i]);
			gaussianBlurring4D.Run();
		}
		cout << "done" << endl;
	}

	// Use identity transformation to start
	mffd = new irtkMultiLevelFreeFormTransformation;
	int tagtrigger = 0;
	// check time sequence is equal
	if(uimage[0]->GetT() != timage[0]->GetT()){
		numberofframeerror();
	}
	for(l=1;l<numberOfUntaggedImages;l++){
		if(uimage[0]->GetT() != uimage[l]->GetT()){
			numberofframeerror();
		}
	}
	for(l=1;l<numberOfTaggedImages;l++){
		if(timage[0]->GetT() != timage[l]->GetT()){
			numberofframeerror();
		}
	}
	for (t = 1; t < uimage[0]->GetT(); t++) {
		tagtrigger = 0;
		// Create registration filter
		irtkCardiac3DImageFreeFormRegistration *cardiacregistration = NULL;
		if (uimage[0]->GetZ() == 1) {
			cerr<<"this mode can't be used with 2D images"<<endl;
			exit(1);
		} 

		// Combine images
		irtkGreyImage **target = new irtkGreyImage*[numberOfUntaggedImages];
		irtkGreyImage **source = new irtkGreyImage*[numberOfUntaggedImages];
		for(l=0;l<numberOfUntaggedImages;l++){
			irtkImageAttributes attr = uimage[l]->GetImageAttributes();
			attr._t = 1;
			target[l] = new irtkGreyImage(attr);
			source[l] = new irtkGreyImage(attr);
			for (z = 0; z < target[l]->GetZ(); z++) {
				for (y = 0; y < target[l]->GetY(); y++) {
					for (x = 0; x < target[l]->GetX(); x++) {
						target[l]->Put(x, y, z, 0, uimage[l]->Get(x, y, z, 0));
						source[l]->Put(x, y, z, 0, uimage[l]->Get(x, y, z, t));
					}
				}
			}
		}

		irtkGreyImage **ttarget = new irtkGreyImage*[numberOfTaggedImages];
		irtkGreyImage **tsource = new irtkGreyImage*[numberOfTaggedImages];
		for(l=0;l<numberOfTaggedImages;l++){
			irtkImageAttributes attr = timage[l]->GetImageAttributes();
			attr._t = 1;
			ttarget[l] = new irtkGreyImage(attr);
			tsource[l] = new irtkGreyImage(attr);
			for (z = 0; z < ttarget[l]->GetZ(); z++) {
				for (y = 0; y < ttarget[l]->GetY(); y++) {
					for (x = 0; x < ttarget[l]->GetX(); x++) {
						ttarget[l]->Put(x, y, z, 0, timage[l]->Get(x, y, z, 0));
						tsource[l]->Put(x, y, z, 0, timage[l]->Get(x, y, z, t));
						//if tag image t and t-1 is identicle use multiple image registration instead.
						if(timage[l]->Get(x,y,z,t) != 0)
							tagtrigger = 1;
					}
				}
			}
		}

        vtkPolyData* tlandmarks = vtkPolyData::New();
        vtkPolyData* slandmarks = vtkPolyData::New();
        if(landmarks != NULL){
            tlandmarks->DeepCopy(landmarks[0]);
            slandmarks->DeepCopy(landmarks[t]);
        }

		if(tagtrigger == 1){
			cardiacregistration = new irtkCardiac3DImageFreeFormRegistration;
			// Set input and output for the registration filter
			cardiacregistration->SetInput(target, source, numberOfUntaggedImages, ttarget, tsource, numberOfTaggedImages);
			if(landmarks != NULL){
				cardiacregistration->SetLandmarks(tlandmarks,slandmarks);
			}
			irtkGreyImage tmpthreshold(*threshold);
			cardiacregistration->SetThreshold(&tmpthreshold);
			cardiacregistration->SetOutput(mffd);
			cardiacregistration->SetDebugFlag(debug);

			// Read parameter if there any, otherwise make an intelligent guess
			cardiacregistration->GuessParameter();
			if (parin_name != NULL) {
				cardiacregistration->irtkMultipleImageRegistration::Read(parin_name);
			}

			// Override parameter settings if necessary
			if (padding != MIN_GREY) {
				cardiacregistration->SetTargetPadding(padding);
			}
			if (spacing > 0) {
				cardiacregistration->SetDX(spacing);
				cardiacregistration->SetDY(spacing);
				cardiacregistration->SetDZ(spacing);
			}

			// Write parameters if necessary
			if (parout_name != NULL) {
				cardiacregistration->irtkMultipleImageRegistration::Write(parout_name);
			}

            if (adaptive > 0) {               
                times = 0;
                t < round(uimage[0]->GetT() / 3) ? 
                    (times = t):(times = round(uimage[0]->GetT() / 3));
                t > round(uimage[0]->GetT()*2 / 3) ? 
                    (times = (uimage[0]->GetT() - 1 - t)):(times = times);
                weight = cardiacregistration->GetLambda1();
                weight = pow(adaptive,2*times)*weight;
                cout << "current lambda1 of frame " << t << " is "<<weight<<endl;
                cardiacregistration->SetLambda1(weight);
                weight = cardiacregistration->GetLambda2();
                weight = pow(adaptive,2*times)*weight;
                cout << "current lambda2 of frame " << t << " is "<<weight<<endl;
                cardiacregistration->SetLambda2(weight);
            }

			// Run registration filter
			cardiacregistration->Run();

			if(cardiacregistration->GetMFFDMode()){
				mffd->CombineLocalTransformation();
			}
		}else{
			cerr << "Empty image?" <<endl;
            exit(1);
		}

		// Write the final transformation estimate
		if (dofout_name != NULL) {
			char buffer[255];
			sprintf(buffer, "%s%d_sequence_%.2d.dof.gz", dofout_name, numberOfTaggedImages, t);
			mffd->irtkTransformation::Write(buffer);
		}
        delete cardiacregistration;

		for (l = 0; l < numberOfUntaggedImages; l++) {
			delete target[l];
			delete source[l];
		}

		delete []target;
		delete []source;

		for (l = 0; l < numberOfTaggedImages; l++) {
			delete ttarget[l];
			delete tsource[l];
		}

		delete []ttarget;
		delete []tsource;
        tlandmarks->Delete();
        slandmarks->Delete();
	}
	if(thresholdname!=NULL)
		delete threshold;
}
