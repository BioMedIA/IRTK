/*=========================================================================

Library   : Image Registration Toolkit (IRTK)
Module    : $Id$
Copyright : Imperial College, Department of Computing
Visual Information Processing (VIP), 2008 onwards
Date      : $Date$
Version   : $Revision$
Changes   : $Author$

Copyright (c) 1999-2014 and onwards, Imperial College London
All rights reserved.
See LICENSE for details

=========================================================================*/

#include <irtkImage.h>
#include <iostream>
#include <irtkTransformation.h>
#include <irtkEMClassification.h>
#include <irtkGaussian.h>

// Default filenames
char *input_name = NULL, *output_name  = NULL;
char *trans_name = NULL;

void usage()
{
    cerr << "Usage: ejecfrac [threshold] [frames] [mask-value] [output] [transformation prefix]\n" << endl;
    cerr << "-landmark [landmark name] landmark information"    << endl;
    exit(1);
}

int main(int argc, char **argv)
{
    int i, j, k, l, t, es, aa, bottom, frames, mask_value, landmarkson, jcbian, ok;
    double x, y, z, minvolume , volume;
    irtkPointSet landmarks;

    landmarkson = 0;

    // Check command line
    if (argc < 6) {
        usage();
    }

    // Parse image
    input_name  = argv[1];
    argc--;
    argv++;
    frames = atoi(argv[1]);
    argc--;
    argv++;
    mask_value = atoi(argv[1]);
    argc--;
    argv++;
    output_name = argv[1];
    argc--;
    argv++;
    trans_name = argv[1];
    argc--;
    argv++;

    while (argc > 1) {
        ok = false;
        if ((ok == false) && (strcmp(argv[1], "-landmark") == 0)) {
            argc--;
            argv++;
            landmarks.ReadVTK(argv[1]);
            landmarkson = 1;
            argc--;
            argv++;
            ok = true;
        }
        if (ok == false) {
            cerr << "Can not parse argument " << argv[1] << endl;
            usage();
        }
    }
    // Read image
    cout << "Reading image ... "; cout.flush();
    irtkRealImage *image = new irtkRealImage(input_name);
    cout << "done" << endl;

    // Preprocess
    if(landmarkson == 1){
        for( i=0; i<landmarks.Size(); i++){
            image->WorldToImage(landmarks(i));
        }
        aa = round(landmarks(5)._z);
        bottom = round(landmarks(6)._z);
        if(aa>bottom) swap(aa,bottom);
        l = 0;
        for( k = aa; k < bottom; k++){
            for( j = 0; j < image->GetY(); j++){
                for( i = 0; i < image->GetX(); i++){
                    if(image->GetAsDouble(i,j,k) == 2 && l == 0){
                        l = k;
                    }
                }
            }
        }
        aa = l;
        if(bottom >= image->GetZ())
            bottom = image->GetZ() - 1;
    }else{
        aa = 0;
        bottom = image->GetZ() - 1;
    }

    minvolume = 0;
    volume = 0;
    int n;
    double m;
    n = 0;

    //irtkTransformation tmp;
    for (k = aa; k < bottom+1; k++) {
        for (j = 0; j < image->GetY(); j++) {
            for (i = 0; i < image->GetX(); i++) {
                if (image->Get(i, j, k) == mask_value) {
                    n ++;
                }
            }
        }
    }
    if( landmarkson == 1)
        volume = n/0.8;
    else
        volume = n;
    minvolume = n;
    es = 1;
    //evaluate min using transformation

    //transform atlas and create tatlas evaluate min using jacobian
    for (t =1; t<frames; t++){
        char buffer[255];
        sprintf(buffer, "%s%.2d.dof.gz", trans_name, t);

        irtkTransformation *mffd1 = NULL;
        irtkTransformation *transform = irtkTransformation::New(buffer);
        if (strcmp(transform->NameOfClass(), "irtkRigidTransformation") == 0) {
            mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkRigidTransformation *)transform));
        } else {
            if (strcmp(transform->NameOfClass(), "irtkAffineTransformation") == 0) {
                mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkAffineTransformation *)transform));
            } else {
                if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
                    mffd1 = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
                } else {
                    cerr << "Input transformation is not of type rigid or affine " << endl;
                    cerr << "or multi-level free form deformation" << endl;
                    exit(1);
                }
            }
        }
        delete transform;

        n = 0;
        m = 0;
        for (k = aa; k < bottom+1; k++) {
            for (j = 0; j < image->GetY(); j++) {
                for (i = 0; i < image->GetX(); i++) {
                    if (image->Get(i, j, k) == mask_value) {
                        x = i;
                        y = j;
                        z = k;
                        image->ImageToWorld(x, y, z);
                        m += mffd1->Jacobian(x, y, z);
                        n ++;
                    }
                }
            }
        }

        if(m<minvolume){
            minvolume = m;
            es = t;
        }

        delete mffd1;
    }

    if( landmarkson == 1)
        minvolume = minvolume/0.8;

    ofstream fout(output_name,ios::app);
    fout << minvolume*image->GetXSize()*image->GetYSize()*image->GetZSize()/1000 << " "
        << volume*image->GetXSize()*image->GetYSize()*image->GetZSize()/1000 << " "
        << (volume - minvolume)/volume << endl << es << endl;
    fout.close();

}
