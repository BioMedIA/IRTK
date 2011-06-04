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

#include <irtkGaussianBlurring.h>

#include <irtkGaussianBlurring4D.h>

char *dofout_name = NULL, *parin_name  = NULL, *parout_name = NULL, *input_name = NULL, *segmentation_name = NULL;

void usage()
{
    cerr << "Usage: motiontrackcardiac2 [image sequence] [segmentation] <options> \n" << endl;
    cerr << "where <options> is one or more of the following:\n" << endl;
    cerr << "<-parin file>        Read parameter from file" << endl;
    cerr << "<-parout file>       Write parameter to file" << endl;
    cerr << "<-dofout file>       Write transformation to file" << endl;
    cerr << "<-ref file>          Reference time frame (default = first frame of image sequence)" << endl;
    cerr << "<-Rx1 value>         Region of interest in images" << endl;
    cerr << "<-Ry1 value>         Region of interest in images" << endl;
    cerr << "<-Rz1 value>         Region of interest in images" << endl;
    cerr << "<-Rt1 value>         Region of interest in images" << endl;
    cerr << "<-Rx2 value>         Region of interest in images" << endl;
    cerr << "<-Ry2 value>         Region of interest in images" << endl;
    cerr << "<-Rz2 value>         Region of interest in images" << endl;
    cerr << "<-Rt2 value>         Region of interest in images" << endl;
    cerr << "<-Tp  value>         Padding value" << endl;
    cerr << "<-lambdac v1...vn>   Lambda coefficients for each segment" << endl;
    exit(1);
}


int main(int argc, char **argv)
{
    int i, n, t, x, y, z, x1, y1, z1, t1, x2, y2, z2, t2, ok, debug;
    double spacing, sigma, xaxis[3], yaxis[3], zaxis[3], *lambdac;
    irtkGreyPixel padding;
    irtkImageFreeFormRegistrationMode mode;
    irtkMultiLevelFreeFormTransformation **mffd;

    // Check command line
    if (argc < 3) {
        usage();
    }

    // Get image names for sequence
    input_name = argv[1];
    argv++;
    argc--;
    segmentation_name = argv[1];
    argv++;
    argc--;

    // Read image sequence
    cout << "Reading image sequence ... "; cout.flush();
    irtkGreyImage *image = new irtkGreyImage(input_name);
    irtkGreyImage *segmentation = new irtkGreyImage(segmentation_name);

    irtkGreyPixel current,min,max,*ptr;
    int sum;
    segmentation->GetMinMax(&min,&max);

    n = 0;
    for (current = min; current <= max; current ++){
        ptr = segmentation->GetPointerToVoxels();
        sum = 0;
        for (i = 0; i < segmentation->GetNumberOfVoxels(); i++){
            if(*ptr == current) sum++;
            ptr++;
        }
        if(sum > 0) n++;
    }

    lambdac = NULL;

    // Fix ROI
    x1 = 0;
    y1 = 0;
    z1 = 0;
    t1 = 0;
    x2 = image->GetX();
    y2 = image->GetY();
    z2 = image->GetZ();
    t2 = image->GetT();

    // Default parameters
    padding   = MIN_GREY;
    spacing   = 0;
    sigma     = 0;
    mode      = RegisterXYZ;
    debug     = false;

    // Parse remaining parameters
    while (argc > 1) {
        ok = false;
        if ((ok == false) && (strcmp(argv[1], "-Rx1") == 0)) {
            argc--;
            argv++;
            x1 = atoi(argv[1]);
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-Rx2") == 0)) {
            argc--;
            argv++;
            x2 = atoi(argv[1]);
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-Ry1") == 0)) {
            argc--;
            argv++;
            y1 = atoi(argv[1]);
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-Ry2") == 0)) {
            argc--;
            argv++;
            y2 = atoi(argv[1]);
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-Rz1") == 0)) {
            argc--;
            argv++;
            z1 = atoi(argv[1]);
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-Rz2") == 0)) {
            argc--;
            argv++;
            z2 = atoi(argv[1]);
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-Rt1") == 0)) {
            argc--;
            argv++;
            t1 = atoi(argv[1]);
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-Rt2") == 0)) {
            argc--;
            argv++;
            t2 = atoi(argv[1]);
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
        if ((ok == false) && (strcmp(argv[1], "-lambdac") == 0)) {
            argc--;
            argv++;
            lambdac = new double[n];
            for(i = 0; i < n; i++){
                lambdac[i] = atof(argv[1]);
                argc--;
                argv++;
            }
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
        if ((ok == false) && (strcmp(argv[1], "-debug") == 0)) {
            argc--;
            argv++;
            ok = true;
            debug = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-xy_only") == 0)) {
            argc--;
            argv++;
            mode = RegisterXY;
            ok = true;
        }

        if (ok == false) {
            cerr << "Can not parse argument " << argv[1] << endl;
            usage();
        }
    }

    // Image orientation
    image->GetOrientation(xaxis, yaxis, zaxis);

    // If there is an region of interest, use it
    if ((x1 != 0) || (x2 != image->GetX()) ||
        (y1 != 0) || (y2 != image->GetY()) ||
        (z1 != 0) || (z2 != image->GetZ()) ||
        (t1 != 0) || (t2 != image->GetT())) {
            *image = image->GetRegion(x1, y1, z1, t1, x2, y2, z2, t2);
    }

    // If sigma is larger than 0, blur images using 4D blurring
    if (sigma > 0) {
        cout << "Blurring image sequences ... "; cout.flush();
        for (i = 0; i < n; i++) {
            irtkGaussianBlurring4D<irtkGreyPixel> gaussianBlurring4D(sigma);
            gaussianBlurring4D.SetInput (image);
            gaussianBlurring4D.SetOutput(image);
            gaussianBlurring4D.Run();
        }
        cout << "done" << endl;
    }

    // Use identity transformation to start
    mffd = new irtkMultiLevelFreeFormTransformation*[n];

    for(i = 0; i < n; i++){
        mffd[i] = new irtkMultiLevelFreeFormTransformation;
    }

    for (t = 1; t < image->GetT(); t++) {

        // Create registration filter
        irtkImageMultiFreeFormRegistration2 *registration = NULL;
        if (image->GetZ() == 1) {
            cerr<<"no 2D for motiontrack2 now"<<endl;
            exit(1);
        } else {
            registration = new irtkImageMultiFreeFormRegistration2;
        }

        // Combine images
        irtkImageAttributes attr = image->GetImageAttributes();
        attr._t = 1;
        irtkGreyImage *target = new irtkGreyImage(attr);
        irtkGreyImage *source = new irtkGreyImage(attr);
        for (z = 0; z < target->GetZ(); z++) {
            for (y = 0; y < target->GetY(); y++) {
                for (x = 0; x < target->GetX(); x++) {
                    target->Put(x, y, z, 0, image->Get(x, y, z, 0));
                    source->Put(x, y, z, 0, image->Get(x, y, z, t));
                }
            }
        }

        // Set input and output for the registration filter
        registration->SetInput(target, source, segmentation);
        registration->SetOutput((irtkTransformation **)mffd);
        registration->SetMode(mode);
        registration->SetDebugFlag(debug);
        if(lambdac != NULL) {
          registration->SetLambdaC(lambdac);
        }

        // Read parameter if there any, otherwise make an intelligent guess
        if (parin_name != NULL) {
            registration->irtkImageRegistration2::Read(parin_name);
        } else {
            registration->GuessParameter();
        }

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
            registration->irtkImageRegistration2::Write(parout_name);
        }

        // Run registration filter
        registration->Run();

        if(registration->GetMFFDMode()){
            for(i = 0; i < n; i++){
            mffd[i]->CombineLocalTransformation();
            }
        }

        // Write the final transformation estimate
        if (dofout_name != NULL) {
            for(i = 0; i < n; i++){
            char buffer[255];
            sprintf(buffer, "%s\\%d_segment_%.2d.dof.gz", dofout_name,i, t);
            mffd[i]->irtkTransformation::Write(buffer);
            }
        }

        delete registration;
        delete target;
        delete source;
    }
    for(i = 0; i < n; i++){
        delete mffd[i];
    }
    delete []mffd;
    if(lambdac != NULL){
        delete []lambdac;
    }
}
