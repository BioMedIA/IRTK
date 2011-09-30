/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2009 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

=========================================================================*/

#include <irtkSegmentationFunction.h>

#include <irtkRegistration2.h>

char *dofout_name = NULL, *dofin_name = NULL, *parin_name  = NULL, *parout_name = NULL, *output_name = NULL, *input_name = NULL;

void usage()
{
    cerr << "Usage: dynamic_segmentation [Input] [Output]" << endl;
    cerr << "<-dofout folder>                       transformation output folder" << endl;
    cerr << "<-atlas n atlas1...atlasn>             atlas information for segmentation" << endl;
    cerr << "<-dofin folder>                        transformation input folder" << endl;
    cerr << "<-parin  file>                         motion tracking parameter file"<< endl;
    cerr << "<-parout file>                         Write parameter to file" << endl;
    cerr << "<-iterations value>                    number of iterations" << endl;
    exit(1);
}

int main( int argc, char** argv )
{
    int i, n, t, x, y, z, ok,numberofatlas;
    irtkMultiLevelFreeFormTransformation *mffd;
    irtkRealImage **atlas = NULL;

    // Check command line
    if (argc < 5) {
        usage();
    }

    n = 1;
    numberofatlas = 0;

    input_name = argv[1];
    argv++;
    argc--;
    output_name = argv[1];
    argv++;
    argc--;

    // Read image sequence
    cout << "Reading image sequence ... "; cout.flush();
    irtkGreyImage *image = new irtkGreyImage(input_name);

    // Parse remaining parameters
    while (argc > 1) {
        ok = false;
        if ((ok == false) && (strcmp(argv[1], "-dofout") == 0)) {
            argc--;
            argv++;
            dofout_name = argv[1];
            argc--;
            argv++;
            ok = true;
        }
        if ((ok == false) && (strcmp(argv[1], "-dofin") == 0)) {
            argc--;
            argv++;
            dofin_name = argv[1];
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
        if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
            argc--;
            argv++;
            ok = true;
            n = atoi(argv[1]);
            argc--;
            argv++;
        }
        if ((ok == false) && (strcmp(argv[1], "-atlas") == 0)) {
            argc--;
            argv++;
            numberofatlas = atoi(argv[1]);
            atlas = new irtkRealImage*[numberofatlas];
            argc--;
            argv++;
            // Read atlas for each tissue
            for (i = 0; i < numberofatlas; i++) {
                atlas[i] = new irtkRealImage;
                atlas[i]->Read(argv[1]);
                cerr << "Image " << i <<" = " << argv[1] <<endl;
                argc--;
                argv++;
            }
            ok = true;
        }
        if (ok == false) {
            cerr << "Can not parse argument " << argv[1] << endl;
            usage();
        }
    }

    if(atlas == NULL){
        cerr << "please set atlas" << endl;
        exit(1);
    }

    irtkImageAttributes attr = image->GetImageAttributes();
    attr._t = 1 + n;
    irtkGreyImage *output = new irtkGreyImage(attr);
    attr._t = 1;
    irtkGreyImage *current = new irtkGreyImage(attr);

    irtkGreyImage *target = new irtkGreyImage(attr);
    for (z = 0; z < target->GetZ(); z++) {
        for (y = 0; y < target->GetY(); y++) {
            for (x = 0; x < target->GetX(); x++) {
                target->PutAsDouble(x, y, z, 0, image->Get(x, y, z, 0));
                current->Put(x,y,z,0);
            }
        }
    }

    // Use identity transformation to start
    mffd = new irtkMultiLevelFreeFormTransformation;

    // first segmentation
    irtkDynamicSegmentation<irtkGreyPixel> dsegmentation;
    dsegmentation.SetInput(target,mffd,atlas,numberofatlas);
    dsegmentation.SetOutput(current);
    dsegmentation.Run();
    for (z = 0; z < current->GetZ(); z++) {
        for (y = 0; y < current->GetY(); y++) {
            for (x = 0; x < current->GetX(); x++) {
                output->PutAsDouble(x,y,z,0,current->GetAsDouble(x,y,z));
            }
        }
    }

    // motiontracking
    for (t = 1; t < image->GetT(); t++) {
        if(dofin_name == NULL){
            // Create registration filter
            irtkImageFreeFormRegistration2 *registration = NULL;
            if (image->GetZ() == 1) {
                cerr<<"no 2D for motiontrack2 now"<<endl;
                exit(1);
            } else {
                registration = new irtkImageFreeFormRegistration2;
            }

            // Combine images
            irtkGreyImage *source = new irtkGreyImage(attr);
            for (z = 0; z < target->GetZ(); z++) {
                for (y = 0; y < target->GetY(); y++) {
                    for (x = 0; x < target->GetX(); x++) {
                        source->PutAsDouble(x, y, z, 0, image->Get(x, y, z, t));
                    }
                }
            }

            // Set input and output for the registration filter
            registration->SetInput(target, source);
            registration->SetOutput(mffd);

      		  // Make an initial Guess for the parameters.
      		  registration->GuessParameter();
      		  // Overrride with any the user has set.
            if (parin_name != NULL) {
                registration->irtkImageRegistration2::Read(parin_name);
            }

            // Write parameters if necessary
            if (parout_name != NULL) {
                registration->irtkImageRegistration2::Write(parout_name);
            }

            // Run registration filter
            registration->Run();

            if(registration->GetMFFDMode()){
                mffd->CombineLocalTransformation();
            }

            // Write the final transformation estimate
            if (dofout_name != NULL) {
                char buffer[255];
                sprintf(buffer, "%ssequence_%.2d.dof.gz", dofout_name, t);
                mffd->irtkTransformation::Write(buffer);
            }

            delete registration;
            delete source;
        }else{
            char buffer[255];
            sprintf(buffer, "%ssequence_%.2d.dof.gz", dofin_name, t);
            cout << "Reading: " << buffer << endl;
            mffd->irtkTransformation::Read(buffer);
        }
        dsegmentation.UpdateMFFD(t);
    }

    for( i = 0; i < n; i ++){

        cout << "Refine segmentation with motion cue iteration: " << i+1 << endl;

        dsegmentation.Run();
        for (z = 0; z < current->GetZ(); z++) {
            for (y = 0; y < current->GetY(); y++) {
                for (x = 0; x < current->GetX(); x++) {
                    output->PutAsDouble(x,y,z,1+i,current->GetAsDouble(x,y,z));
                }
            }
        }
    }       
    output->Write(output_name);
    delete mffd;
    delete target;
    delete image;
    delete output;
    if(atlas != NULL){
         for (i = 0; i < numberofatlas; i++) {
           delete atlas[i];
         }
         delete []atlas;
    }
    delete current;
}
