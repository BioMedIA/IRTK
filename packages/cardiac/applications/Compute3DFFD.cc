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

#include <irtkRegistration2.h>

void usage()
{
    cerr << "Usage: compute3DFFD [dofin] [dofoutPrefix] [N] <options> \n" << endl;
    cerr << "takes in a 4D transformation and computes the [N] 3D transformation for time points n*1/N \n" << endl;
    cerr << "<-invert>             Invert the transformation" << endl;

    exit(1);
}

int main(int argc, char **argv)
{
    double t, dt;
    int N, n, Nx, Ny, Nz, i, ok, invert;
    double x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis[3], yaxis[3], zaxis[3];
    double *cpX = NULL, *cpY = NULL, *cpZ = NULL, 
        *dispX = NULL, *dispY = NULL, *dispZ = NULL,
        *tmpdispX = NULL, *tmpdispY = NULL, *tmpdispZ = NULL,
        error, errorOld;

    int maxIter = 100;

    // Default filenames
    char *dofin_name  = NULL, *dofout_prefix = NULL, *dofout_current = NULL;

    invert = false;

    // Check command line
    if (argc < 4) {
        usage();
    }

    // Parse source and target images
    dofin_name = argv[1];
    argc--;
    argv++;
    dofout_prefix = argv[1];
    argc--;
    argv++;
    N = atoi(argv[1]);
    argc--;
    argv++;

    while (argc > 1) {
        ok = false;
        if ((ok == false) && (strcmp(argv[1], "-invert") == 0)) {
            argc--;
            argv++;
            invert = true;
            ok = true;
        }
        if (ok == false) {
            cerr << "Can not parse argument " << argv[1] << endl;
            usage();
        }
    }

    dt = 1./double(N);
    irtkBSplineFreeFormTransformation *ffd = NULL;

    // Read target image
    cout << "Reading transformation ... "; cout.flush();
    irtkTransformation *transform = irtkTransformation::New(dofin_name);
    irtkMultiLevelFreeFormTransformation *mffd = NULL;
    irtkBSplineFreeFormTransformationPeriodic *tffd = NULL;
    if (strcmp(transform->NameOfClass(), "irtkMultiLevelFreeFormTransformation") == 0) {
        mffd = new irtkMultiLevelFreeFormTransformation(*((irtkMultiLevelFreeFormTransformation *)transform));
    } else {
        cerr << "Transformation has to be of type irtkMultiLevelFreeFormTransformation" << endl;
        exit(1);
    }
    // get finest level
    tffd = (irtkBSplineFreeFormTransformationPeriodic *)
        mffd->GetLocalTransformation(mffd->NumberOfLevels() - 1);
    cout << "done" << endl;

    // initialisation parameters
    Nx = tffd->GetX();
    Ny = tffd->GetY();
    Nz = tffd->GetZ();
    x1 = y1 = z1 = 0;
    tffd->LatticeToWorld(x1, y1, z1);
    x2 = Nx-1;
    y2 = Ny-1;
    z2 = Nz-1;
    tffd->LatticeToWorld(x2, y2, z2);
    dx = tffd->GetXSpacing()/2;
    dy = tffd->GetYSpacing()/2;
    dz = tffd->GetZSpacing()/2;
    tffd->GetOrientation(xaxis, yaxis, zaxis);

    Nx = Nx*2;
    Ny = Ny*2;
    Nz = Nz*2;

    // control points
    cpX = new double[Nx*Ny*Nz];
    cpY = new double[Nx*Ny*Nz];
    cpZ = new double[Nx*Ny*Nz];

    cout << "computing 3D transformations ... " << endl; cout.flush();
    for (n = 0; n < N; n++) {
        cout<<n;
        t = double(n)*dt;

        // initialise 3D FFD
        ffd = new irtkBSplineFreeFormTransformation(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis);

        error = 10;

        dispX = new double[Nx*Ny*Nz];
        dispY = new double[Nx*Ny*Nz];
        dispZ = new double[Nx*Ny*Nz];

        tmpdispX = new double[Nx*Ny*Nz];
        tmpdispY = new double[Nx*Ny*Nz];
        tmpdispZ = new double[Nx*Ny*Nz];

        //create control point locations;
        for (i = 0; i < Nx*Ny*Nz; i++) {
            ffd->ControlPointLocation(i, cpX[i], cpY[i], cpZ[i]);
            dispX[i] = cpX[i];
            dispY[i] = cpY[i];
            dispZ[i] = cpZ[i];
            mffd->LocalDisplacement(dispX[i], dispY[i], dispZ[i], t);
        }
        if (invert == true) {
            for (i = 0; i < Nx*Ny*Nz; i++) {      
                cpX[i] += dispX[i];
                cpY[i] += dispY[i];
                cpZ[i] += dispZ[i];

                dispX[i] = -dispX[i];
                dispY[i] = -dispY[i];
                dispZ[i] = -dispZ[i];
            }
        }

        //////////////////////////////////////////
        // run Approximate function multiple times
        int iter = 0;
        do {

            for (i = 0; i < Nx*Ny*Nz; i++) {
                tmpdispX[i] = dispX[i];
                tmpdispY[i] = dispY[i];
                tmpdispZ[i] = dispZ[i];
            }

            ffd->Approximate(cpX, cpY, cpZ, tmpdispX, tmpdispY, tmpdispZ, Nx*Ny*Nz);

            errorOld = error;

            // calculate error
            error = 0;
            for (i = 0; i < Nx*Ny*Nz; i++) {
                error += sqrt( tmpdispX[i]*tmpdispX[i] + tmpdispY[i]*tmpdispY[i] + tmpdispZ[i]*tmpdispZ[i]);
            }
            error /= 3*Nx*Ny*Nz;

            //cout<<" e: "<<error<<endl;
        } while(errorOld - error > 0.000001 && iter < maxIter); // check for convergence
        //////////////////////////////////////////
        //////////////////////////////////////////
        delete []tmpdispX;
        delete []tmpdispY;
        delete []tmpdispZ;
        delete []dispX;
        delete []dispY;
        delete []dispZ;

        cout<<" -> final error: "<<error<<endl;

        // write new ffd at time t
        irtkMultiLevelFreeFormTransformation *mtffd = new irtkMultiLevelFreeFormTransformation;
        mtffd->PushLocalTransformation(ffd);
        dofout_current = new char[200];
        sprintf(dofout_current, "%s_t%d.dof.gz", dofout_prefix, n);
        mtffd->irtkTransformation::Write(dofout_current);

        delete dofout_current;
        delete mtffd;
    }
    cout << "done" << endl;

    delete mffd;
}
